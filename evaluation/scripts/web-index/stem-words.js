#!env node

var readline = require('readline');
var natural = require('natural');

var rl = readline.createInterface({
  input: process.stdin,
  output: process.stdout,
  terminal: false
});

rl.on('line', function (line) {
  //console.log(natural.PorterStemmer.stem(line));
  console.log(natural.LancasterStemmer.stem(line));
});
