use aws_sdk_dynamodb::config::Region;
use aws_sdk_dynamodb::types::AttributeValue::S;
use aws_sdk_dynamodb::types::{
    AttributeDefinition, BillingMode, KeySchemaElement, KeyType, ScalarAttributeType,
};
use aws_sdk_dynamodb::Client;
use std::collections::{HashMap, HashSet};

const TABLENAME: &str = "rdv";

pub async fn create_rdv_table_if_not_exists(client: &Client) {
    let cmd = client
        .list_tables()
        .send()
        .await
        .expect("Failed to list tables");
    if !cmd.table_names.unwrap().contains(&TABLENAME.to_owned()) {
        create_rdv_table(client).await;
    }
}

pub async fn create_rdv_table(client: &Client) {
    let ad = AttributeDefinition::builder()
        .attribute_name("key")
        .attribute_type(ScalarAttributeType::S)
        .build()
        .unwrap();

    let ks = KeySchemaElement::builder()
        .attribute_name("key")
        .key_type(KeyType::Hash)
        .build()
        .unwrap();

    let cmd = client
        .create_table()
        .table_name(TABLENAME)
        .key_schema(ks)
        .attribute_definitions(ad)
        .billing_mode(BillingMode::PayPerRequest);
    cmd.send().await.unwrap();
}

pub async fn delete_rdv_table(client: &Client) {
    let cmd = client.delete_table().table_name(TABLENAME);
    cmd.send().await.unwrap();
}

pub async fn delete_key(client: &Client, key: &str) {
    let cmd = client
        .delete_item()
        .table_name(TABLENAME)
        .key("key", S(key.to_owned()));
    cmd.send().await.unwrap();
}

pub async fn add_attr(client: &Client, key: &str, attr_key: &str, attr_value: &str) {
    let update = client
        .update_item()
        .table_name(TABLENAME)
        .key("key", S(key.to_owned()))
        .expression_attribute_names("#attr_key", attr_key)
        .expression_attribute_values(":attr_value", S(attr_value.to_owned()))
        .update_expression("SET #attr_key = :attr_value");
    update.send().await.unwrap();
}

pub async fn get_attr(client: &Client, key: &str, attr_key: &str) -> Option<String> {
    let cmd = client
        .get_item()
        .table_name(TABLENAME)
        .key("key", S(key.to_owned()))
        .projection_expression(attr_key)
        .consistent_read(true);
    let res = cmd.send().await.unwrap();
    if res.item.is_none() {
        return None;
    }
    match res.item.unwrap().get(attr_key) {
        Some(S(s)) => Some(s.clone()),
        _ => None,
    }
}

pub async fn add_attr_if_not_exist(
    client: &Client,
    key: &str,
    attr_key: &str,
    attr_value: &str,
) -> String {
    let update = client
        .update_item()
        .table_name(TABLENAME)
        .key("key", S(key.to_owned()))
        .expression_attribute_names("#attr_key", attr_key)
        .expression_attribute_values(":attr_value", S(attr_value.to_owned()))
        .condition_expression("attribute_not_exists(#attr_key)")
        .update_expression("SET #attr_key = :attr_value");
    let res = update.send().await;
    let update_success = res.is_ok();
    if update_success {
        attr_value.to_owned()
    } else {
        let cmd = client
            .get_item()
            .table_name(TABLENAME)
            .key("key", S(key.to_owned()))
            .projection_expression(attr_key)
            .consistent_read(true);
        let res = cmd.send().await.unwrap();
        match res.item.unwrap().get(attr_key) {
            Some(S(s)) => s.clone(),
            _ => panic!("must exist"),
        }
    }
}

// attr_keys is a mutable reference to a set of strings.
// if the key is found, it is removed from the set.
pub async fn get_attrs(
    client: &Client,
    key: &str,
    attr_keys: &mut HashSet<String>,
) -> HashMap<String, String> {
    let proj = attr_keys
        .iter()
        .map(|k| k.to_owned())
        .collect::<Vec<_>>()
        .join(", ");
    let cmd = client
        .get_item()
        .table_name(TABLENAME)
        .key("key", S(key.to_owned()))
        .projection_expression(proj)
        .consistent_read(true);
    let res = cmd.send().await.unwrap();
    let mut vs: HashMap<String, String> = HashMap::new();
    for (k, v) in res.item.unwrap().iter() {
        attr_keys.remove(k);
        vs.insert(k.to_owned(), v.as_s().unwrap().to_owned());
    }
    vs
}

pub async fn make_db_client() -> Client {
    let config = aws_config::defaults(aws_config::BehaviorVersion::latest())
        .region(Region::new("us-east-1"))
        .load()
        .await;
    Client::new(&config)
}
