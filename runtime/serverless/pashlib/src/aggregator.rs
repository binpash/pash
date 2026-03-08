use anyhow::{anyhow, Result};
use std::io::{copy, stdout, Write};
use std::process::{Command, Stdio};

pub trait Aggregator {
    fn aggregate(&self, inputs: &[String], output: &str) -> Result<()>;
}

pub struct SortMergeAggregator;

impl Aggregator for SortMergeAggregator {
    fn aggregate(&self, inputs: &[String], output: &str) -> Result<()> {
        if inputs.is_empty() {
            return Ok(());
        }
        if inputs.len() == 1 {
            if output == "-" {
                let mut src = std::fs::File::open(&inputs[0])?;
                let mut out = stdout();
                copy(&mut src, &mut out)?;
                out.flush()?;
            } else {
                std::fs::copy(&inputs[0], output)?;
            }
            return Ok(());
        }

        let mut cmd = Command::new("sort");
        cmd.arg("-m");
        for p in inputs {
            cmd.arg(p);
        }
        cmd.stdout(Stdio::piped());

        let mut child = cmd.spawn()?;
        let mut child_out = child
            .stdout
            .take()
            .ok_or_else(|| anyhow!("failed to capture sort -m stdout"))?;

        if output == "-" {
            let mut out = stdout();
            copy(&mut child_out, &mut out)?;
            out.flush()?;
        } else {
            let mut out = std::fs::File::create(output)?;
            copy(&mut child_out, &mut out)?;
            out.flush()?;
        }

        let status = child.wait()?;
        if !status.success() {
            return Err(anyhow!("sort -m failed with status {}", status));
        }

        Ok(())
    }
}
