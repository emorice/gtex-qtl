"""
VCF input utilities
"""

import struct
import gzip

def _read_exact(stream, length):
    """
    Read exactly `length` bytes from `stream` by making several calls to `read` if necessary
    """
    buf = stream.read(length)
    if len(buf) < length:
        if not buf:
            raise EOFError
        # This should be very rarely needed
        return buf + read_exact(stream, length - len(buf))
    return buf

# BGZFixedHeader = namedtuple('GZFixedHeader', 'id1 id2 compression_method flags '
#        'mtime extra_flags operating_system extra_length extra_id1 extra_id2 '
#        'extra_length_inner block_size')
_BGZFFormat = 'BBBBIBB' 'H' 'BBH' 'H'
_BGZFSize = struct.calcsize(_BGZFFormat)

class _GZFlags:
    FEXTRA = 1 << 2

def _read_bgzip_block_size(stream):
    """
    Read BGZIP header and extract block size. This will consume some bytes
    """
    (
        id1, id2, _cm, flags, _mtime, _xflags, _os,
        _xlen,
        extra_id1, extra_id2, extra_length_inner, block_size
        ) = struct.unpack(_BGZFFormat, _read_exact(stream, _BGZFSize))

    if id1 != 31 or id2 != 139: # GZIP magic
        raise ValueError('Not a GZIP file')
    if not flags & _GZFlags.FEXTRA:
        raise ValueError('Not a BGZIP file: missing extra fields')
    # BGZIP magic with two bytes of payload
    if (extra_id1, extra_id2, extra_length_inner) != (66, 67, 2):
        # The spec allows a BGZIP file to contain other extra fields inserted
        # before the BGZIP extra field proper. To be implemented the day we
        # encounter such a file.
        raise NotImplementedError('Unknown extra field')
    return block_size + 1

def _read_next_coo(stream):
    """
    Extract next coordinate, when stream is at the start of a valid gzip block
    """
    with gzip.open(stream, 'rt', encoding='ascii') as ustream:
        # Discard likely truncated line
        ustream.readline()
        # Discard any comments
        while True:
            line = ustream.readline()
            if line[0] != '#':
                break

        return line.split('\t', maxsplit=2)[:2]

def simple_vcf_index(file, chunk_size=4<<20):
    with open(file, 'rb') as stream:
        block_size = 0
        pos = 0
        next_checkpoint_pos = 0
        index = []
        while block_size != 28: # Size of EOF marker block
            block_size = _read_bgzip_block_size(stream)
            if pos >= next_checkpoint_pos:
                stream.seek(pos)
                index.append((
                        *_read_next_coo(stream), pos
                        ))
                next_checkpoint_pos += chunk_size
            stream.seek(pos + block_size)
            pos += block_size
        return index
