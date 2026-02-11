use pash_sls_lib::db_helper::*;
use pash_sls_lib::holepunch::*;
use tokio::io;
use tokio::fs::File;
use tokio::net::TcpStream;
use tokio::io::AsyncReadExt;
use tokio::io::AsyncWriteExt;

pub async fn setup(me: &str, peer: &str, rdv_key: &str) -> TcpStream {
    let client = make_db_client().await;
    create_rdv_table_if_not_exists(&client).await;
    let mut ctx = PashCtx::new(me, rdv_key).await;
    let stream = ctx.connect(peer).await;
    stream
}

#[tokio::main]
async fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut handles = Vec::new();
    for arg in args[1..].iter() {
        if arg.starts_with("send") {
            let v: Vec<&str> = arg.split("*").collect();
            let rdv_key = v[1].to_string();
            let me = v[2].to_string();
            let peer = v[3].to_string();
            let fifo_name = v[4].to_string();
            let handle = tokio::spawn(async move {
                // println!("send: me={}, peer={}, rdv_key={}, fifo_name={}", me, peer, rdv_key, fifo_name);
                let stream = setup(&me, &peer, &rdv_key).await;
                let (mut rd, mut wr) = stream.into_split();
                let mut buf = [0u8; 128];
                let mut fifo = File::open(fifo_name).await.unwrap();
                rd.read_exact(&mut buf).await.unwrap();
                // copy from fifo to wr
                io::copy(&mut fifo, &mut wr).await.unwrap();
            });
            handles.push(handle);
        } else {
            let v: Vec<&str> = arg.split("*").collect();
            let rdv_key = v[1].to_string();
            let me = v[2].to_string();
            let peer = v[3].to_string();
            let fifo_name = v[4].to_string();
            let handle = tokio::spawn(async move {
                // println!("recv: me={}, peer={}, rdv_key={}, fifo_name={}", me, peer, rdv_key, fifo_name);
                let stream = setup(&me, &peer, &rdv_key).await;
                let (mut rd, mut wr) = stream.into_split();
                let buf = [0u8; 128];
                let mut fifo = File::create(fifo_name).await.unwrap();
                wr.write_all(&buf).await.unwrap();
                io::copy(&mut rd, &mut fifo).await.unwrap();
            });
            handles.push(handle);
        }
    }

    for handle in handles {
        handle.await.unwrap();
    }
}
