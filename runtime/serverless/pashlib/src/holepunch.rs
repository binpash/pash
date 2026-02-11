use serde::{Deserialize, Serialize};
use tokio::net::{TcpSocket, TcpStream};
use crate::db_helper::*;
use crate::stun_helper;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PashCtx {
    #[serde(skip)]
    pub rdv_key: String,
    #[serde(skip)]
    pub name: String,
    pub local_addr: String,
    pub external_addr: String,
}

impl PashCtx {
    pub async fn new(name: &str, key: &str) -> Self {
        Self {
            rdv_key: key.to_owned(),
            name: name.to_owned(),
            local_addr: "".to_owned(),
            external_addr: "".to_owned(),
        }
    }

    pub fn attr_key(name: &str) -> String {
        "N_".to_owned() + name
    }

    pub async fn connect(&mut self, dst_name: &str) -> TcpStream {
        // decide ec2 vs lambda based on env variable
        if std::env::var("LEASH_ON_LAMBDA").is_ok() {
            self.connect_lambda(dst_name).await
        } else {
            self.connect_ec2(dst_name).await
        }
    }

    pub async fn connect_ec2(&mut self, dst_name: &str) -> TcpStream {
        let client = make_db_client().await;

        // 0) perform stun to get external addr
        let (local_addr, external_addr) = stun_helper::get_addr().await;
        self.local_addr = local_addr.to_owned();
        self.external_addr = external_addr.to_owned();

        // 1) create tcp socket using an ephemeral port
        let sock = TcpSocket::new_v4().unwrap();
        sock.bind("0.0.0.0:0".parse().unwrap()).unwrap(); // 0 = let OS pick
        let local_addr = sock.local_addr().unwrap();
        let ip = self.external_addr.parse::<std::net::SocketAddr>().unwrap().ip();
        let share_addr = format!(
            "{}:{}",
            ip.to_string(),
            local_addr.port().to_string()
        );

        // 2) share my info
        let info_str = serde_json::to_string(&Self {
            rdv_key: "".to_owned(),
            name: "".to_owned(),
            local_addr: self.local_addr.clone(),
            external_addr: share_addr,
        })
        .unwrap();
        add_attr_if_not_exist(&client, &self.rdv_key, &Self::attr_key(&self.name), &info_str).await;

        // 3) check peer info
        let info_str = loop {
            let res = get_attr(&client, &self.rdv_key, &Self::attr_key(dst_name)).await;
            if let Some(info_str) = res {
                break info_str;
            }
            tokio::time::sleep(tokio::time::Duration::from_secs(1)).await;
        };
        let dst_info: PashCtx = serde_json::from_str(&info_str).unwrap();

        // 4) connect to peer
        let stream = sock
            .connect(dst_info.external_addr.parse().unwrap())
            .await
            .unwrap();
        stream.set_nodelay(true).unwrap();
        stream
    }

    pub async fn connect_lambda(&mut self, dst_name: &str) -> TcpStream {
        let client = make_db_client().await;

        // 0) perform stun to get external addr; if binding failed, retry
        const MAX_ATTEMPTS: usize = 10;
        const RETRY_DELAY_MS: u64 = 10;

        let mut attempt = 0usize;
        let mut last_err: Option<anyhow::Error>;

        let (sock, local_addr, external_addr) = loop {
            attempt += 1;

            let res: anyhow::Result<(TcpSocket, String, String)> = async {
                let (local_addr_global, external_addr_global) = stun_helper::get_addr().await;
                let local_port = local_addr_global
                    .parse::<std::net::SocketAddr>()
                    .unwrap()
                    .port();

                let sock = TcpSocket::new_v4().unwrap();
                sock.set_reuseport(true).unwrap();

                let bind_addr: std::net::SocketAddr = format!("0.0.0.0:{}", local_port)
                    .parse()
                    .unwrap();
                sock.bind(bind_addr).unwrap();

                Ok((sock, bind_addr.to_string(), external_addr_global))
            }
            .await;

            match res {
                Ok(v) => break v,
                Err(e) => {
                    last_err = Some(e);
                    if attempt >= MAX_ATTEMPTS {
                        panic!("Stun/bind failed after {} attempts: {:?}", attempt, last_err);
                    }
                    tokio::time::sleep(tokio::time::Duration::from_millis(RETRY_DELAY_MS)).await;
                }
            }
        };

        // 1) share my info
        let info_str = serde_json::to_string(&Self {
            rdv_key: "".to_owned(),
            name: "".to_owned(),
            local_addr: local_addr.clone(),
            external_addr: external_addr.clone(),
        })
        .unwrap();
        add_attr_if_not_exist(&client, &self.rdv_key, &Self::attr_key(&self.name), &info_str).await;

        // 2) check peer info
        let info_str = loop {
            let res = get_attr(&client, &self.rdv_key, &Self::attr_key(dst_name)).await;
            if let Some(info_str) = res {
                break info_str;
            }
            tokio::time::sleep(tokio::time::Duration::from_secs(1)).await;
        };
        let dst_info: PashCtx = serde_json::from_str(&info_str).unwrap();

        // 3) connect to peer
        let stream = sock
            .connect(dst_info.external_addr.parse().unwrap())
            .await
            .unwrap();
        stream.set_nodelay(true).unwrap();
        stream
    }
}
