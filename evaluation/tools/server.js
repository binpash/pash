#!env node
const http = require('http');
const fs = require('fs');
const cp = require('child_process');
const crypto = require('crypto');
const HASH = "f4f8b9dc33fe805bcd0a854b19554f1daa596cf46ae3ef811a27b7d490dae42c";
const PORT = 6001

let jobs = [];

let validateObj = (obj) => {
  let keys = ["dt", "program"];
  for (var i = 0; i < keys.length; i++) {
    if (!obj[keys[i]]) {
      return {success: false, msg: "object missing " + keys[i]}
    }
  };

  let h = crypto.createHash('sha256').update(obj.dt).digest('hex');
  if (h !== HASH) {
    return {success: false, msg: "security token failed"}
  }
  return {success: true};
}

let validateURI = (uri) => {
  const parsed = parseInt(uri, 10);
  if (isNaN(parsed)) {
    return {success: false, msg: "Did not call `/<int>`"};
  } 
  if ((parsed < 1) || (parsed > jobs.length - 1)) {
    return {success: false, msg: "No such ID"};
  }
  return {success: true, msg: data[parsed]};
}

let newJob = (obj) => {
  console.log(jobs, jobs.length, obj);
  if (!validateObj(obj).success) {
    return validateObj(obj)
  }
  const jobID = jobs.length;
  jobs[jobID] = obj;
  const filename = jobID + ".sh";
  try {
    fs.writeFileSync(filename, obj.program, "utf-8");
  } catch (e) {
    return e
  }
//  const pipeline = `
//  nc -l ${obj.workerPort} |
//        ./${fileName} |
//         nc -C ${serverIP} ${serverPort}`
  let c = cp.exec(`sh ./${filename}`, (error, stdout, stderr) => {
    if (error) {
      console.error(`exec error: ${error}`);
      jobs[jobID].error = error;
      return
//    cb({success: false});
    }
    jobs[jobID].complete = true;
    if (stdout) {
      jobs[jobID].stdout = stdout;
    }
    if (stderr) {
      jobs[jobID].stderr = stderr;
    }
  });
  // return a second after the process started
  return {success: true}
}

http.createServer( (request, response) => {
  if (request.method === 'PUT' && request.url === '/') {
    // PUT /control adds to the global timeseries
    let body = [];
    request.on('data', (chunk) => {
      body.push(chunk);
    }).on('end', () => {
      body = Buffer.concat(body).toString();
      console.log(body);
      response.end(JSON.stringify(newJob(JSON.parse(body))));
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
