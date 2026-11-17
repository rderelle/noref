//! Run log: the lines printed during a run are also written to `<prefix>.log`, without
//! the progress bars, after the command line and followed by the run time or by the
//! error that stopped the run.
use anyhow::Context;
use std::fmt::Display;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;
use std::time::Instant;

/// Print a line (formatted like `eprintln!`) and write it to the run log.
macro_rules! say {
    ($log:expr, $($arg:tt)*) => {
        $log.line(format_args!($($arg)*))
    };
}
pub(crate) use say;

struct LogFile {
    writer: BufWriter<File>,
    /// First write error, reported when the log is closed.
    error: Option<io::Error>,
}

pub(crate) struct RunLog {
    file: Mutex<LogFile>,
    started: Instant,
}

impl RunLog {
    /// Create the log file and write the command line to it.
    pub(crate) fn create(path: &Path) -> anyhow::Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("cannot create the log file {}", path.display()))?;
        let log = Self {
            file: Mutex::new(LogFile {
                writer: BufWriter::new(file),
                error: None,
            }),
            started: Instant::now(),
        };
        let command: Vec<String> = std::env::args_os()
            .map(|arg| arg.to_string_lossy().into_owned())
            .collect();
        log.write(format_args!("command: {}", command.join(" ")));
        Ok(log)
    }

    /// Print a line to the standard error and write it to the log (see [`say`]).
    pub(crate) fn line(&self, text: impl Display) {
        eprintln!("{text}");
        self.write(&text);
    }

    /// Write a line to the log only.
    fn write(&self, text: impl Display) {
        let mut file = self.file.lock().expect("log mutex poisoned");
        if file.error.is_none() {
            if let Err(e) = writeln!(file.writer, "{text}") {
                file.error = Some(e);
            }
        }
    }

    /// Close the log after the run time, or after the error that stopped the run.
    pub(crate) fn finish(self, result: &anyhow::Result<()>) -> anyhow::Result<()> {
        match result {
            Ok(()) => self.write(format_args!(
                "total time: {:.1}s",
                self.started.elapsed().as_secs_f64()
            )),
            Err(e) => self.write(format_args!("Error: {e:#}")),
        }
        let mut file = self.file.into_inner().expect("log mutex poisoned");
        if let Some(e) = file.error {
            return Err(e).context("cannot write the log file");
        }
        file.writer.flush().context("cannot write the log file")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lines_are_logged_after_the_command_and_before_the_outcome() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.log");
        let log = RunLog::create(&path).unwrap();
        say!(log, "# step {}", 1);
        say!(log, " . count: {}", 42);
        log.finish(&Ok(())).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert!(lines[0].starts_with("command: "), "{text}");
        assert_eq!(&lines[1..3], ["# step 1", " . count: 42"]);
        assert!(
            lines[3].starts_with("total time: ") && lines.len() == 4,
            "{text}"
        );

        let log = RunLog::create(&path).unwrap();
        say!(log, "# step 1");
        log.finish(&Err(anyhow::anyhow!("no input"))).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.ends_with("# step 1\nError: no input\n"), "{text}");
    }

    #[test]
    fn a_log_in_a_missing_directory_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let error = RunLog::create(&dir.path().join("missing/run.log"))
            .err()
            .unwrap();
        assert!(error.to_string().contains("cannot create the log file"));
    }
}
