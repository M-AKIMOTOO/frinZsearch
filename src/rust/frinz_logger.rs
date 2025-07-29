use std::fs::{File, OpenOptions};
use std::io::{self, Write};
use std::sync::{Arc, Mutex};

pub struct Logger {
    console_output_enabled: bool,
    file_output_enabled: bool,
    file_handle: Option<Arc<Mutex<File>>>,
}

impl Logger {
    pub fn new() -> Self {
        Logger {
            console_output_enabled: true, // Console output enabled by default
            file_output_enabled: false,
            file_handle: None,
        }
    }

    pub fn setup(&mut self, console_enabled: bool, file_enabled: bool, filepath: Option<&str>) -> io::Result<()> {
        self.console_output_enabled = console_enabled;
        self.file_output_enabled = file_enabled;

        if self.file_output_enabled {
            if let Some(path) = filepath {
                let file = OpenOptions::new()
                    .create(true)
                    .write(true)
                    .append(true) // Append mode
                    .open(path)?;
                self.file_handle = Some(Arc::new(Mutex::new(file)));
            } else {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, "File output is enabled, but no file path was specified."));
            }
        }
        Ok(())
    }

    pub fn log(&self, message: &str) {
        if self.console_output_enabled {
            println!("{}", message);
        }
        if self.file_output_enabled {
            if let Some(file_handle) = &self.file_handle {
                let mut file = file_handle.lock().unwrap();
                if let Err(e) = writeln!(file, "{}", message) {
                    eprintln!("Error writing to file: {}", e);
                }
            }
        }
    }

    // Helper for logging formatted strings
    pub fn log_fmt(&self, args: std::fmt::Arguments) {
        self.log(&format!("{}", args));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Read;
    use tempfile::tempdir;

    #[test]
    fn test_logger_console_output() {
        // This test cannot directly capture console output,
        // so it only verifies that the log function runs without error
        let logger = Logger::new();
        logger.log("Test message to console");
        logger.log_fmt(format_args!("Formatted test message: {}", 123));
    }

    #[test]
    fn test_logger_file_output() -> io::Result<()> {
        let dir = tempdir()?;
        let file_path = dir.path().join("test_log.txt");
        let file_path_str = file_path.to_str().unwrap();

        let mut logger = Logger::new();
        logger.setup(false, true, Some(file_path_str))?;

        logger.log("First log entry.");
        logger.log_fmt(format_args!("Second log entry: {}", "data"));

        // ファイルの内容を読み込んで検証
        let mut file = File::open(&file_path)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;

        assert!(contents.contains("First log entry."));
        assert!(contents.contains("Second log entry: data"));
        assert_eq!(contents.lines().count(), 2);

        Ok(())
    }

    #[test]
    fn test_logger_no_output() -> io::Result<()> {
        let dir = tempdir()?;
        let file_path = dir.path().join("test_log_no_output.txt");
        let file_path_str = file_path.to_str().unwrap();

        let mut logger = Logger::new();
        // Disable both console and file output
        logger.setup(false, false, Some(file_path_str))?;

        logger.log("This should not appear anywhere.");

        // Verify that the file does not exist or is empty
        assert!(!file_path.exists() || fs::read_to_string(&file_path)?.is_empty());

        Ok(())
    }
}
