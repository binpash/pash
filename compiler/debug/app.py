from flask import Flask, render_template
from werkzeug.exceptions import abort

app = Flask(__name__)
worker_db = {}

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/<int:worker_id>')
def get_worker_log(worker_id):
    return render_template('wlog.html')

def log():
    addr = request.remote_addr
    return 


class Worker:
    def __init__(self, id):
        self.id = id
        self.log = ""

    def clear_log(self):
        self.log = ""