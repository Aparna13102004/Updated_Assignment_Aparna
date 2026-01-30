import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from plasmid_builder import build_plasmid


def run():
    if len(sys.argv) < 3:
        sys.stderr.write("Expected arguments: <input_fasta> <design_file>\n")
        sys.exit(1)

    fasta_path, design_path = sys.argv[1:3]

    assembled_sequence = build_plasmid(fasta_path, design_path)

    plasmid_record = SeqRecord(
        seq=Seq(assembled_sequence),
        id="Output_Plasmid",
        description="Plasmid constructed using host ORI and design file"
    )

    output_path = "Output.fa"
    SeqIO.write(plasmid_record, output_path, format="fasta")

    print(f"Plasmid sequence successfully saved to {output_path}")


if __name__ == "__main__":
    run()

