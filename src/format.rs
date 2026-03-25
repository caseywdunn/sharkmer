/// Format a count with comma-separated thousands (e.g. 1,234,567).
pub(crate) fn format_count(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result
}

/// Format bytes as a human-readable string (e.g. "1.2 GB").
pub(crate) fn format_bytes(bytes: u64) -> String {
    const KB: f64 = 1024.0;
    const MB: f64 = KB * 1024.0;
    const GB: f64 = MB * 1024.0;
    let b = bytes as f64;
    if b < KB {
        format!("{} B", bytes)
    } else if b < MB {
        format!("{:.1} KB", b / KB)
    } else if b < GB {
        format!("{:.1} MB", b / MB)
    } else {
        format!("{:.1} GB", b / GB)
    }
}

/// Format a duration as a human-readable string (e.g. "1m 23s").
pub(crate) fn format_duration(d: std::time::Duration) -> String {
    let secs = d.as_secs_f64();
    if secs < 60.0 {
        format!("{:.1}s", secs)
    } else if secs < 3600.0 {
        let mins = (secs / 60.0).floor() as u64;
        let remaining = secs - (mins as f64 * 60.0);
        format!("{}m {:.0}s", mins, remaining)
    } else {
        let hours = (secs / 3600.0).floor() as u64;
        let remaining = secs - (hours as f64 * 3600.0);
        let mins = (remaining / 60.0).floor() as u64;
        format!("{}h {}m", hours, mins)
    }
}
