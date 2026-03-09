use pash_sls_lib::recovery::FtExecutor;
use tracing::Level;

fn log_level_from_env() -> Level {
    let v = std::env::var("RUST_LOG").unwrap_or_else(|_| "info".to_string());
    match v.to_ascii_lowercase().as_str() {
        "trace" => Level::TRACE,
        "debug" => Level::DEBUG,
        "warn" => Level::WARN,
        "error" => Level::ERROR,
        _ => Level::WARN,
    }
}

#[tokio::main]
async fn main() {
    tracing_subscriber::fmt()
        .with_max_level(log_level_from_env())
        .with_target(false)
        .init();

    let args: Vec<String> = std::env::args().collect();
    let mut handles = Vec::new();

    for arg in args[1..].iter() {
        let executor = FtExecutor::new(arg).unwrap_or_else(|e| {
            panic!("failed to parse endpoint arg '{}': {}", arg, e);
        });

        let handle = tokio::spawn(async move {
            let _ = executor.execute().await.unwrap();
        });
        handles.push(handle);
    }

    for handle in handles {
        handle.await.unwrap();
    }
}
