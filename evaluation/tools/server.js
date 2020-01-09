#!env node
const http = require('http');
const fs = require('fs');
const cp = require('child_process');
const crypto = require('crypto');
const HASH = "f4f8b9dc33fe805bcd0a854b19554f1daa596cf46ae3ef811a27b7d490dae42c";
const PORT = 6001

// nc -l 5000 | tr '[:lower:]' '[:upper:]' | nc -C 158.130.4.212 5555

let jobs = [];

// example object:
// {dt: "passphrase", workerPort: 1234, serverIP: "192.168.0.1", serverPort: 4321, program: "..."  }
let validateObj = (obj) => {
  let keys = ["dt", "workerPort", "serverIP", "serverPort", "program"];
  for (var i = 0; i < keys.length; i++) {
    if (!obj[keys[i]]) {
      return {"success": false, "msg": "object missing " + keys[i]}
    }
  };

  let h = crypto.createHash('sha256').update(obj.dt).digest('hex');
  if (h !== HASH) {
    return {"success": false, "msg": "security token failed"}
  }
}

let validateURI = (uri) => {
  const parsed = parseInt(uri, 10);
  if (isNaN(parsed)) {
    return {"success": false, "msg": "Did not call `/<int>`"};
  } 
  if ((parsed < 1) || (parsed > jobs.lenth - 1)) {
    return {"success": false, "msg": "No such ID"};
  }
  reutn {"success": true, "msg": data[parsed]};
}

let newJob = (obj) => {
  const jobID = jobs.lenth;
  const fileName = jobID + ".sh";
  fs.writeFileSync(fileName, obj.program, "utf-8");
  const pipeline = `
  nc -l ${obj.workerPort} |
        ./${fileName} |
         nc -C ${serverIP} ${serverPort}`
  let c = exec(pipeline, (error, stdout, stderr) => {
    if (error) {
      console.error(`exec error: ${error}`);
      return;
    }
    console.log(`stdout: ${stdout}`);
    console.error(`stderr: ${stderr}`);
  });
}

function roughScale(x, base) {
  const parsed = parseInt(x, base);
  if (isNaN(parsed)) { return 0 }
  return parsed * 100;
}

http.createServer( (request, response) => {
  if (request.method === 'PUT' && request.url === '/newjob') {
    // PUT /control adds to the global timeseries
    let body = [];
    request.on('data', (chunk) => {
      body.push(chunk);
    }).on('end', () => {
      body = Buffer.concat(body).toString();
      response.end(JSON.stringify(newJob(JSON.parse(body)));
      );
    });
  } else if (request.url === '/') {
    // otherwise, return all current jobs
    response.writeHead(200);
    response.end(JSON.stringify(jobs));
  } else {
    let v = validateURI(request.url);
		response.writeHead(v.success? 200 : 404);
		response.end(JSON.stringify(v));
  }
}).listen(6000);
