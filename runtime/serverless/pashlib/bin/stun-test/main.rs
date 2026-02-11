use pash_sls_lib::stun_helper;
use tokio::io;

#[tokio::main]
async fn main() {
    let (local_addr, external_addr) = stun_helper::get_addr().await;
    println!("local_addr: {}, external_addr: {}", local_addr, external_addr);
}
