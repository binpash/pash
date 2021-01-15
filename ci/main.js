#! /usr/bin/env node

const { createGitHubWebhookServer } = require('./github-webhooks.js');

if (require.main === module) {
    createGitHubWebhookServer(2047);
}
