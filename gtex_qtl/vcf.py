"""
VCF input utilities
"""

import struct
import logging
import gzip
import numpy as np

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

        chrom, pos = line.split('\t', maxsplit=2)[:2]
        return chrom, int(pos)

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

class VCFFields:
    """
    Fixed fields for genotype files according to VCF spec
    """
    CHROM = 0
    POS = 1
    ID = 2
    REF = 3
    ALT = 4
    QUAL = 5
    FILTER = 6
    INFO = 7
    FORMAT = 8
    SAMPLES = 9

def _dosage_from_records(records):
    if not records:
        return np.empty((0, 0), dtype=bool), np.empty((0, 0), dtype=np.uint8)
    # Compact everything as an aligned matrix of bytes
    # Note that this will detect any line with a wrong number of bytes
    gt_data = np.stack([np.frombuffer(rec[VCFFields.SAMPLES].encode('ascii'), dtype='u1') for rec in records])
    # Discard field separators
    gt_data = gt_data[:, ::2]

    # Reshape
    n_sites, n_alleles = gt_data.shape
    assert not n_alleles % 2
    n_samples = n_alleles // 2
    gt_data = gt_data.reshape((n_sites, n_samples, 2))

    # Missing mask
    missing = gt_data == np.uint8(ord('.'))
    missing = missing[..., 0] | missing[..., 1]

    # Make numerical
    gt_data -= np.uint8(ord('0'))

    # Sum
    gt_data = np.add(gt_data[..., 0], gt_data[..., 1], dtype=np.uint8)

    # At this point, missing entries will contain unspecified other values. We
    # leave it as such since it makes easier to detect errors

    # Checks
    assert np.all(missing | (gt_data <= 2))

    return missing, gt_data

def _read_region(file, chrom, start, end, skip=0):
    """
    Args:
        start, end: closed interval
    """
    lines = []
    in_lines = iter(gzip.open(file, 'rt', encoding='utf8'))
    for _ in range(skip):
        next(in_lines)
    chrom_found = False
    for line in in_lines:
        if line[0] == '#':
            continue
        # 8 Fixed fields, plus a mandatory FORMAT field
        fields = line.split('\t', maxsplit=9)
        if fields[VCFFields.CHROM] != chrom:
            if chrom_found:
                return
            continue
        chrom_found = True
        pos = int(fields[VCFFields.POS])
        fields[VCFFields.POS] = pos
        if pos < start:
            continue
        if pos > end:
            return
        yield fields

def _parse_region(file, chrom, start, end, skip=0):
    records = list(_read_region(file, chrom, start, end, skip=0))
    meta = [rec[:VCFFields.SAMPLES] for rec in records]
    missing, dosage = _dosage_from_records(records)
    return meta, missing, dosage

def parse_region_indexed(file, index, chrom, start, end, skip=0):
    """
    Faster variant of parse_region when a simple index is available to seek
    close to the wanted region
    """
    chrom_found = False
    last_before = index[0]
    for entry in index:
        i_chrom, i_start, i_pos = entry
        if i_chrom != chrom:
            if chrom_found:
                break
            last_before = entry
        else:
            chrom_found = True
            if i_start >= start:
                break
            last_before = entry
    if not chrom_found:
        logging.warning('Could not find %s in index, falling back to unindexed',
                chrom)
        last_before = index[0]

    file_pos = last_before[-1]
    # If starting from the middle of the file, discard one truncated line
    skip = 0 if file_pos == 0 else 1

    stream = open(file, 'rb')
    stream.seek(file_pos)
    return _parse_region(stream, chrom, start, end, skip=skip)
