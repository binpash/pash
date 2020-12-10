#!env node
const http = require('http');
const fs = require('fs');

const workers = {
  "local": {ip: "127.0.0.1", port: 6000},
  "beta":  {ip: "158.130.4.113", port: 6000},
  "gamma": {ip: "158.130.4.114", port: 6000},
  "delta": {ip: "158.130.4.120", port: 6000}
};

// ./client.js            -- reports global status
const help=`
./client.js <w>        -- reports on status of worker <w>
./client <w> <expr>    -- launches <expr> on <w>

Workers currently available:
  ${Object.keys(workers).map( (e) => e ).join('\n  ')}
`;

let get = (key, cb) => {
  http.get(`http://${workers[key].ip}:${workers[key].port}/`, (res) => {
    res.setEncoding('utf8');
    let rawData = '';
    res.on('data', (chunk) => { rawData += chunk; });
    res.on('end', () => {
      cb(rawData);
    });
  }).on('error', (e) => {
    cb(`Got error: ${e.message}`);
  });
};

let generateHeader = (key, data, op) => {
  return {
    hostname: workers[key].ip,
    port: workers[key].port,
    path: '/',
    method: op || 'PUT',
    headers: {
      'Content-Type': 'application/json',
      'Content-Length': data.length
    }
  }
}

let put = (key, json, cb) => {
  const req = http.request(generateHeader(key, json), (res) => {
    console.log(`statusCode: ${res.statusCode}`)

    res.on('data', (d) => {
      process.stdout.write(d)
    })
  })
  
  req.on('error', (error) => {
    cb(error)
  })
  
  req.end(json)
  //req.end() // TODO combine
}
  
const dtVar="PASH_TOKEN";

let printAndExit = (msg) => {
  console.log(msg);
  process.exit(-1);
}

let checkDishToken = () => {
  if (!process.env[dtVar]) {
    printAndExit(`Need to set ${dtVar} environment variable`);
  }
  return process.env[dtVar]; 
}

let checkArgs = (args) => {
  if (args.length === 3 && /help/.test(args[2])) {
    printAndExit(help);
  }
  if (args.length === 3 && !workers[args[2]]) {
    printAndExit(`${args[2]} does not exist in workers!`);
  }
}

checkArgs(process.argv)
if (process.argv.length < 3) {
    printAndExit(help);
} else if (process.argv.length === 3) {
  let w = workers[process.argv[2]]
  get(process.argv[2], console.log)
} else {
  let program;
  if (process.argv[3] === "-f") {
    program = fs.readFileSync(process.argv[4], 'utf-8');
  } else {
    program = process.argv[3]
  }
  if (program.indexOf('\'')) {
    //FIXME Confirm it does not matter
    console.error('WARNING, program includes single quotes that might mess things up');
  }
  let job = JSON.stringify({ dt: checkDishToken(), program: program })
  console.log(job)
  put(process.argv[2], job, console.log);
}
// checkDishToken()
// if (process.argv[1])
// if (process.arb
// get(/*construct url*/, console.log)
