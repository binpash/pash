use std::sync::Arc;

use stun::agent::TransactionId;
use stun::client::ClientBuilder;
use stun::message::{Getter, Message, BINDING_REQUEST, BINDING_SUCCESS};
use stun::xoraddr::XorMappedAddress;
use tokio::net::UdpSocket;

pub async fn get_addr() -> (String, String) {
    let server = "stun.l.google.com:19302";

    let (handler_tx, mut handler_rx) = tokio::sync::mpsc::unbounded_channel();

    let conn = UdpSocket::bind("0:0").await.unwrap();
    conn.connect(server).await.unwrap();

    let local_addr = conn.local_addr().unwrap().to_string();

    let cc = Arc::from(conn);
    let mut client = ClientBuilder::new().with_conn(cc).build().unwrap();

    let mut msg = Message::new();
    let txn_id = TransactionId::new();
    msg.build(&[Box::new(txn_id), Box::new(BINDING_REQUEST)])
        .unwrap();

    client.send(&msg, Some(Arc::new(handler_tx))).await.unwrap();

    let event = handler_rx.recv().await.unwrap();
    let msg = event.event_body.unwrap();
    assert_eq!(msg.typ, BINDING_SUCCESS);
    assert_eq!(msg.transaction_id, txn_id);
    let mut xor_addr = XorMappedAddress::default();
    xor_addr.get_from(&msg).unwrap();

    let addr = xor_addr.to_string();

    client.close().await.unwrap();

    (local_addr, addr)
}
