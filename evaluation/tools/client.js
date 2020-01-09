#!env node

const workers = [
  {name: "local", ip: "127.0.0.1", port: 6000},
  {name: "beta",  ip: "158.130.4.113", port: 6000},
  {name: "gamma", ip: "158.130.4.114", port: 6000},
  {name: "delta", ip: "158.130.4.120", port: 6000}
];

// ./client.js            -- reports global status
const help=`
./client.js <w>        -- reports on status of worker <w>
./client <w> <expr>    -- launches <expr> on <w>

Workers currently available:
  ${workers.map( (e) => e.name ).join('\n  ')}
`;

let get = (url, cb) => {
  http.get(url, (res) => {
    res.setEncoding('utf8');
    let rawData = '';
    res.on('data', (chunk) => { rawData += chunk; });
    res.on('end', () => {
      cb(rawData);
    });
  }).on('error', (e) => {
    cb(`Got error: ${e.message}`);
  });
}

const dtVar="DISH_TOKEN";

let printAndExit = (msg) => {
  console.log(msg);
  process.exit(-1);
}

let checkDishToken = () => {
  if (!process.env[dtVar]) {
    printAndExit(`Need to set ${dtVar} environment variable`);
  }
}

let checkArgs = (args) => {
  if (args.length == 3 && /help/.test(args[2])) {
    printAndExit(help);
  }
  if (args.length == 3 && !workers[args[2]]) {
    printAndExit(`${args[2]} does not exist in workers!`);
  }
}

checkArgs(process.argv)
if (process.argv.lenth < 3) {
    printAndExit(help);
} else if (process.argv.lenth === 3) {
  let w = workers[process.argv[2]]
  get(`http://${w.ip}:${w.port}/`, console.log)
} else {
  checkDishToken();
  console.error(`implement PUT`);
}
// checkDishToken()
// if (process.argv[1])
// if (process.arb
// get(/*construct url*/, console.log)
