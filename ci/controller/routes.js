// A controller implementation that uses AWS

const { respond, getRequestBody } = require('./lib.js');
const SSM = require('@aws-sdk/client-ssm');

const ssm = new SSM.SSMClient({ region: "us-east-1" });

const echo = async (req, res) =>
      respond(res, 200, await getRequestBody(req));

const assignment = (assignees) =>
    new SSM.SendCommandCommand({
        DocumentName: 'AWS-RunShellScript',
        DocumentVersion: '1',
        Targets: [{Key: 'tag:Pash', Values: [assignees] }],
        OutputS3BucketName: 'pash-reports',
        Parameters: {
            commands: ['/home/ubuntu/worker-script.sh'],
            workingDirectory: [''],
            executionTimeout: ['3600']
        },
    });


// Tells every EC2 instance responsible for correctness tests to run their scripts.
const ci = async (req, res) => {
    try {
        const { '$metadata': { httpStatusCode }, Command } = await ssm.send(assignment('ci-correctness'));

        if (httpStatusCode === 200) {
            const { CommandId }  = Command;
            respond(res, 200, `CI signal sent. SSM Command Id: ${CommandId}\n`);
        } else {
            respond(res, 500, `AWS SSM responded with error code ${httpStatusCode}\n`);
        }
    } catch (e) {
        err(e);
        respond(res, 500, `Failed to distribute work to ci workers.\n`);
    }
};


// Reports instance history
// TODO: Upgrade to HTML response so that links can take user to exact command pages in AWS.
const now = async (req, res) => {
    try {
        const command = new SSM.ListCommandInvocationsCommand({ details: true });

        const { '$metadata': { httpStatusCode }, CommandInvocations: history } = await ssm.send(command);

        if (httpStatusCode === 200) {
            res.writeHead(200, { 'Content-Type': 'text/plain' });

            const groupedHistory = history.reduce((res, item) => {
                const { InstanceName: id } = item;
                const existing = res[id] || [];
                return Object.assign(res, { [id]: existing.concat([item]) });
            }, {});

            const sortedHistory = Object.entries(groupedHistory).sort(([a,],[b,]) => b.localeCompare(a));

            for (const [id, instanceHistory] of sortedHistory) {
                if (id === '') continue;

                res.write(`Instance ${id}\n`);
                for (const { CommandId, RequestedDateTime, Status } of instanceHistory) {
                    res.write(`${RequestedDateTime.toISOString()} ${Status} ${CommandId}\n`);
                }
                res.write(`\n\n`);
            }

            res.end();
        } else {
            respond(res, 500, `AWS SSM responded with error code ${httpStatusCode}`);
        }
    } catch (e) {
        console.error(e);
    }
};


module.exports = {
    '/now': now,
    '/echo': echo,
    '/ci': ci,
};
