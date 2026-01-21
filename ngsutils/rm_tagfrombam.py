import pysam
import os

path_in = os.path.join(os.environ.get('HOME'), "tmp/test.bam")
path_out = "./out.bam"

alignments_in = pysam.AlignmentFile(path_in, "rb")
alignments_out = pysam.AlignmentFile(path_out, "wb", template=alignments_in)

for (i, r) in enumerate(alignments_in):
    r.cigarstring = "150M"
    r.set_tag("NM", None)
    alignments_out.write(r)
    if i > 5:
        break

alignments_out.close()
alignments_out.close()
