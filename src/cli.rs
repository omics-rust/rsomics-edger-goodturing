use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Section};

use rsomics_edger_goodturing::{Opts, run};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-edger-goodturing", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    pub counts: PathBuf,
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    #[arg(long)]
    proportions: bool,
    #[arg(long, default_value_t = 1.96)]
    conf: f64,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        let opts = Opts {
            proportions: self.proportions,
            conf: self.conf,
        };
        let n = run(&self.counts, &opts, &mut out)?;
        if !self.common.quiet {
            eprintln!("{n} genes Good-Turing smoothed across columns");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Simple Good-Turing frequency smoothing of count vectors into proportions (edgeR goodTuring/goodTuringProportions).",
    origin: None,
    usage_lines: &["<counts.tsv> [--proportions] [--conf 1.96] [-o out.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "proportions",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Emit a same-shape proportion matrix (goodTuringProportions); zeros share the unseen mass and each column sums to 1.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "conf",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("1.96"),
                description: "Confidence multiplier for the Turing→linear-Good-Turing switch rule.",
                why_default: Some("edgeR's default conf."),
            },
        ],
    }],
    examples: &[
        Example {
            description: "Per-distinct-count smoothed proportions + P0 per column",
            command: "rsomics-edger-goodturing counts.tsv -o gt.tsv",
        },
        Example {
            description: "Same-shape proportion matrix",
            command: "rsomics-edger-goodturing counts.tsv --proportions -o prop.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
