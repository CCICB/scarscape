use std::fmt;

/// Trinucleotide Mutation Counts
///
/// Describes every possible trinculoetide mutation (pyrimidine-centered)
#[derive(Default, Debug)]
pub struct CountsTriChange {
    aca_aaa: u64,
    acc_aac: u64,
    acg_aag: u64,
    act_aat: u64,
    aca_aga: u64,
    acc_agc: u64,
    acg_agg: u64,
    act_agt: u64,
    aca_ata: u64,
    acc_atc: u64,
    acg_atg: u64,
    act_att: u64,
    ata_aaa: u64,
    atc_aac: u64,
    atg_aag: u64,
    att_aat: u64,
    ata_aca: u64,
    atc_acc: u64,
    atg_acg: u64,
    att_act: u64,
    ata_aga: u64,
    atc_agc: u64,
    atg_agg: u64,
    att_agt: u64,
    cca_caa: u64,
    ccc_cac: u64,
    ccg_cag: u64,
    cct_cat: u64,
    cca_cga: u64,
    ccc_cgc: u64,
    ccg_cgg: u64,
    cct_cgt: u64,
    cca_cta: u64,
    ccc_ctc: u64,
    ccg_ctg: u64,
    cct_ctt: u64,
    cta_caa: u64,
    ctc_cac: u64,
    ctg_cag: u64,
    ctt_cat: u64,
    cta_cca: u64,
    ctc_ccc: u64,
    ctg_ccg: u64,
    ctt_cct: u64,
    cta_cga: u64,
    ctc_cgc: u64,
    ctg_cgg: u64,
    ctt_cgt: u64,
    gca_gaa: u64,
    gcc_gac: u64,
    gcg_gag: u64,
    gct_gat: u64,
    gca_gga: u64,
    gcc_ggc: u64,
    gcg_ggg: u64,
    gct_ggt: u64,
    gca_gta: u64,
    gcc_gtc: u64,
    gcg_gtg: u64,
    gct_gtt: u64,
    gta_gaa: u64,
    gtc_gac: u64,
    gtg_gag: u64,
    gtt_gat: u64,
    gta_gca: u64,
    gtc_gcc: u64,
    gtg_gcg: u64,
    gtt_gct: u64,
    gta_gga: u64,
    gtc_ggc: u64,
    gtg_ggg: u64,
    gtt_ggt: u64,
    tca_taa: u64,
    tcc_tac: u64,
    tcg_tag: u64,
    tct_tat: u64,
    tca_tga: u64,
    tcc_tgc: u64,
    tcg_tgg: u64,
    tct_tgt: u64,
    tca_tta: u64,
    tcc_ttc: u64,
    tcg_ttg: u64,
    tct_ttt: u64,
    tta_taa: u64,
    ttc_tac: u64,
    ttg_tag: u64,
    ttt_tat: u64,
    tta_tca: u64,
    ttc_tcc: u64,
    ttg_tcg: u64,
    ttt_tct: u64,
    tta_tga: u64,
    ttc_tgc: u64,
    ttg_tgg: u64,
    ttt_tgt: u64,
    other: u64,
}

/// Fetch Middle character from a string slice
fn middle_char(s: &str) -> Option<char> {
    let mut chars = s.chars();
    let mid = chars.clone().count() / 2;
    chars.nth(mid)
}

fn reverse_complement(seq: &str) -> &str {
    todo!("write a reverse complement script")
}

/// Expects sequences to be lowercase. We should fix that
fn pyrimidine_center(reference: &str, alternative: &str) -> (&str, &str) {
    let nchars = reference.chars().count();

    if nchars % 2 != 0 {
        log::warn!("Can NOT pyrimidine-center a sequence with an even length ({reference})");
        return (&reference, &alternative);
    }

    // Grab middle character from reference sequence
    let Some(middlebase) = middle_char(reference) else {
        // If middle_char returns None return reference and alternative seqs as is
        return (&reference, &alternative);
    };

    // If middle char is a purine, reverse complement both ref and alt to get mutation in terms of
    // pyrimidine ref
    if ['a', 'g'].contains(&middlebase) {
        // reverse complement both sequences
        todo!("Add reverse complement to reference and alternative and return")
    }

    // Otherwise return as is
    (&reference, &alternative)
}

impl CountsTriChange {
    /// Increment the counter for this 3-mer.
    /// Converts trinucleotide str to lowercase, pyrimidine-centered forms (T,C)
    /// by reverse complementing purine-centered sequences.
    /// Will count non-ATCG containing sequences as 'other'
    pub fn increment(&mut self, tri_ref: &str, tri_alt: &str) {
        // Change to lowercase
        let tri_ref_lower = tri_ref.to_ascii_lowercase();
        let tri_alt_lower = tri_alt.to_ascii_lowercase();

        // Change ref to pyridine centered
        // middle base of reference sequence must be C or T;
        // if it is not then reverse complement both ref and alt seqs to ensure it is.

        match (tri_ref_lower.as_str(), tri_alt_lower.as_str()) {
            // Pyrimidine
            ("aca", "aaa") => self.aca_aaa += 1,
            ("acc", "aac") => self.acc_aac += 1,
            ("acg", "aag") => self.acg_aag += 1,
            ("act", "aat") => self.act_aat += 1,
            ("aca", "aga") => self.aca_aga += 1,
            ("acc", "agc") => self.acc_agc += 1,
            ("acg", "agg") => self.acg_agg += 1,
            ("act", "agt") => self.act_agt += 1,
            ("aca", "ata") => self.aca_ata += 1,
            ("acc", "atc") => self.acc_atc += 1,
            ("acg", "atg") => self.acg_atg += 1,
            ("act", "att") => self.act_att += 1,
            ("ata", "aaa") => self.ata_aaa += 1,
            ("atc", "aac") => self.atc_aac += 1,
            ("atg", "aag") => self.atg_aag += 1,
            ("att", "aat") => self.att_aat += 1,
            ("ata", "aca") => self.ata_aca += 1,
            ("atc", "acc") => self.atc_acc += 1,
            ("atg", "acg") => self.atg_acg += 1,
            ("att", "act") => self.att_act += 1,
            ("ata", "aga") => self.ata_aga += 1,
            ("atc", "agc") => self.atc_agc += 1,
            ("atg", "agg") => self.atg_agg += 1,
            ("att", "agt") => self.att_agt += 1,
            ("cca", "caa") => self.cca_caa += 1,
            ("ccc", "cac") => self.ccc_cac += 1,
            ("ccg", "cag") => self.ccg_cag += 1,
            ("cct", "cat") => self.cct_cat += 1,
            ("cca", "cga") => self.cca_cga += 1,
            ("ccc", "cgc") => self.ccc_cgc += 1,
            ("ccg", "cgg") => self.ccg_cgg += 1,
            ("cct", "cgt") => self.cct_cgt += 1,
            ("cca", "cta") => self.cca_cta += 1,
            ("ccc", "ctc") => self.ccc_ctc += 1,
            ("ccg", "ctg") => self.ccg_ctg += 1,
            ("cct", "ctt") => self.cct_ctt += 1,
            ("cta", "caa") => self.cta_caa += 1,
            ("ctc", "cac") => self.ctc_cac += 1,
            ("ctg", "cag") => self.ctg_cag += 1,
            ("ctt", "cat") => self.ctt_cat += 1,
            ("cta", "cca") => self.cta_cca += 1,
            ("ctc", "ccc") => self.ctc_ccc += 1,
            ("ctg", "ccg") => self.ctg_ccg += 1,
            ("ctt", "cct") => self.ctt_cct += 1,
            ("cta", "cga") => self.cta_cga += 1,
            ("ctc", "cgc") => self.ctc_cgc += 1,
            ("ctg", "cgg") => self.ctg_cgg += 1,
            ("ctt", "cgt") => self.ctt_cgt += 1,
            ("gca", "gaa") => self.gca_gaa += 1,
            ("gcc", "gac") => self.gcc_gac += 1,
            ("gcg", "gag") => self.gcg_gag += 1,
            ("gct", "gat") => self.gct_gat += 1,
            ("gca", "gga") => self.gca_gga += 1,
            ("gcc", "ggc") => self.gcc_ggc += 1,
            ("gcg", "ggg") => self.gcg_ggg += 1,
            ("gct", "ggt") => self.gct_ggt += 1,
            ("gca", "gta") => self.gca_gta += 1,
            ("gcc", "gtc") => self.gcc_gtc += 1,
            ("gcg", "gtg") => self.gcg_gtg += 1,
            ("gct", "gtt") => self.gct_gtt += 1,
            ("gta", "gaa") => self.gta_gaa += 1,
            ("gtc", "gac") => self.gtc_gac += 1,
            ("gtg", "gag") => self.gtg_gag += 1,
            ("gtt", "gat") => self.gtt_gat += 1,
            ("gta", "gca") => self.gta_gca += 1,
            ("gtc", "gcc") => self.gtc_gcc += 1,
            ("gtg", "gcg") => self.gtg_gcg += 1,
            ("gtt", "gct") => self.gtt_gct += 1,
            ("gta", "gga") => self.gta_gga += 1,
            ("gtc", "ggc") => self.gtc_ggc += 1,
            ("gtg", "ggg") => self.gtg_ggg += 1,
            ("gtt", "ggt") => self.gtt_ggt += 1,
            ("tca", "taa") => self.tca_taa += 1,
            ("tcc", "tac") => self.tcc_tac += 1,
            ("tcg", "tag") => self.tcg_tag += 1,
            ("tct", "tat") => self.tct_tat += 1,
            ("tca", "tga") => self.tca_tga += 1,
            ("tcc", "tgc") => self.tcc_tgc += 1,
            ("tcg", "tgg") => self.tcg_tgg += 1,
            ("tct", "tgt") => self.tct_tgt += 1,
            ("tca", "tta") => self.tca_tta += 1,
            ("tcc", "ttc") => self.tcc_ttc += 1,
            ("tcg", "ttg") => self.tcg_ttg += 1,
            ("tct", "ttt") => self.tct_ttt += 1,
            ("tta", "taa") => self.tta_taa += 1,
            ("ttc", "tac") => self.ttc_tac += 1,
            ("ttg", "tag") => self.ttg_tag += 1,
            ("ttt", "tat") => self.ttt_tat += 1,
            ("tta", "tca") => self.tta_tca += 1,
            ("ttc", "tcc") => self.ttc_tcc += 1,
            ("ttg", "tcg") => self.ttg_tcg += 1,
            ("ttt", "tct") => self.ttt_tct += 1,
            ("tta", "tga") => self.tta_tga += 1,
            ("ttc", "tgc") => self.ttc_tgc += 1,
            ("ttg", "tgg") => self.ttg_tgg += 1,
            ("ttt", "tgt") => self.ttt_tgt += 1,
            _ => self.other += 1,
        }
    }

    /// Total mutations
    pub fn total(&self, include_other: bool) -> u64 {
        let mut total = self.aca
            + self.acc
            + self.acg
            + self.act
            + self.ata
            + self.atc
            + self.atg
            + self.att
            + self.cca
            + self.ccc
            + self.ccg
            + self.cct
            + self.cta
            + self.ctc
            + self.ctg
            + self.ctt
            + self.gca
            + self.gcc
            + self.gcg
            + self.gct
            + self.gta
            + self.gtc
            + self.gtg
            + self.gtt
            + self.tca
            + self.tcc
            + self.tcg
            + self.tct
            + self.tta
            + self.ttc
            + self.ttg
            + self.ttt;

        if include_other {
            total += self.other;
        }

        total
    }

    /// Core printer: writes a two-column table with the given delimiter
    pub fn fmt_with_delimiter(&self, f: &mut fmt::Formatter<'_>, delim: char) -> fmt::Result {
        // header
        writeln!(f, "context{d}count", d = delim)?;

        // contexts in the same order as your fields
        let pairs = [
            ("ACA", self.aca),
            ("ACC", self.acc),
            ("ACG", self.acg),
            ("ACT", self.act),
            ("ATA", self.ata),
            ("ATC", self.atc),
            ("ATG", self.atg),
            ("ATT", self.att),
            ("CCA", self.cca),
            ("CCC", self.ccc),
            ("CCG", self.ccg),
            ("CCT", self.cct),
            ("CTA", self.cta),
            ("CTC", self.ctc),
            ("CTG", self.ctg),
            ("CTT", self.ctt),
            ("GCA", self.gca),
            ("GCC", self.gcc),
            ("GCG", self.gcg),
            ("GCT", self.gct),
            ("GTA", self.gta),
            ("GTC", self.gtc),
            ("GTG", self.gtg),
            ("GTT", self.gtt),
            ("TCA", self.tca),
            ("TCC", self.tcc),
            ("TCG", self.tcg),
            ("TCT", self.tct),
            ("TTA", self.tta),
            ("TTC", self.ttc),
            ("TTG", self.ttg),
            ("TTT", self.ttt),
            ("other", self.other),
        ];

        for &(ctx, cnt) in &pairs {
            writeln!(f, "{ctx}{d}{cnt}", d = delim)?;
        }
        Ok(())
    }
}

impl fmt::Display for CountsTri {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // default to tab-delimited
        self.fmt_with_delimiter(f, '\t')
    }
}
