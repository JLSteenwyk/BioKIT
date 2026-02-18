import textwrap

from .base import Genome


class GenomeAssemblyMetrics(Genome):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        output_format = self.normalize_output_format(self.output_format)
        l50 = self.calc_l50()
        l90 = self.calc_l90()
        n50 = self.calc_n50()
        n90 = self.calc_n90()

        char_freq = self._get_base_counts()
        assembly_size = self.sum_of_scaffold_lengths()

        if assembly_size == 0:
            gc_content = 0.0
            a_freq = 0.0
            t_freq = 0.0
            g_freq = 0.0
            c_freq = 0.0
        else:
            gc_content = (char_freq.get("G", 0) + char_freq.get("C", 0)) / assembly_size
            a_freq = char_freq.get("A", 0) / assembly_size
            t_freq = char_freq.get("T", 0) / assembly_size
            g_freq = char_freq.get("G", 0) / assembly_size
            c_freq = char_freq.get("C", 0) / assembly_size

        num_of_scaffolds = self.number_of_scaffolds()

        longest_scaffold = self.longest_scaffold()

        (
            num_of_large_scaffolds,
            len_of_large_scaffolds,
        ) = self.number_of_large_scaffolds()

        summary = {
            "Assembly size": assembly_size,
            "L50": l50,
            "L90": l90,
            "N50": n50,
            "N90": n90,
            "GC content": round(gc_content, 4),
            "Number of scaffolds": num_of_scaffolds,
            "Number of large scaffolds": num_of_large_scaffolds,
            "Sum length of large scaffolds": len_of_large_scaffolds,
            "Longest scaffold": longest_scaffold,
            "Frequency of A": round(a_freq, 4),
            "Frequency of T": round(t_freq, 4),
            "Frequency of C": round(c_freq, 4),
            "Frequency of G": round(g_freq, 4),
        }

        if output_format == "tsv":
            print(
                textwrap.dedent(
                    f"\
                {summary['Assembly size']}\tAssembly size\n\
                {summary['L50']}\tL50\n\
                {summary['L90']}\tL90\n\
                {summary['N50']}\tN50\n\
                {summary['N90']}\tN90\n\
                {summary['GC content']}\tGC content\n\
                {summary['Number of scaffolds']}\tNumber of scaffolds\n\
                {summary['Number of large scaffolds']}\tNumber of large scaffolds\n\
                {summary['Sum length of large scaffolds']}\tSum length of large scaffolds\n\
                {summary['Longest scaffold']}\tLongest scaffold\n\
                {summary['Frequency of A']}\tFrequency of A\n\
                {summary['Frequency of T']}\tFrequency of T\n\
                {summary['Frequency of C']}\tFrequency of C\n\
                {summary['Frequency of G']}\tFrequency of G"
                )
            )
            return

        structured_summary = {
            "assembly_size": summary["Assembly size"],
            "l50": summary["L50"],
            "l90": summary["L90"],
            "n50": summary["N50"],
            "n90": summary["N90"],
            "gc_content": summary["GC content"],
            "number_of_scaffolds": summary["Number of scaffolds"],
            "number_of_large_scaffolds": summary["Number of large scaffolds"],
            "sum_length_of_large_scaffolds": summary["Sum length of large scaffolds"],
            "longest_scaffold": summary["Longest scaffold"],
            "frequency_of_a": summary["Frequency of A"],
            "frequency_of_t": summary["Frequency of T"],
            "frequency_of_c": summary["Frequency of C"],
            "frequency_of_g": summary["Frequency of G"],
        }
        print(self.format_object(structured_summary, output_format))

    def process_args(self, args) -> dict:
        if args.threshold is None:
            threshold = 500
        else:
            threshold = int(args.threshold)

        return dict(
            fasta=args.fasta,
            threshold=threshold,
            output_format=getattr(args, "format", None),
        )
