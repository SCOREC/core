# How to use
# 1. Login to https://app.globus.org/settings/developers and copy a project app id and secret
# 2. Use the id and secret to create and endpoint https://funcx.readthedocs.io/en/latest/sdk.html#client-credentials-with-clients
#     $ export FUNCX_SDK_CLIENT_ID="b0500dab-ebd4-430f-b962-0c85bd43bdbb"
#     $ export FUNCX_SDK_CLIENT_SECRET="ABCDEFGHIJKLMNOP0123456789="
# 3. Set up an endpoint on the computer that will run the tests, using these instructions: https://funcx.readthedocs.io/en/latest/endpoints.html
# 4. Create install-test.sh and run-test.sh on target computer

from globus_compute_sdk import Executor
import sys
import os

name = sys.argv[1]
branch = sys.argv[2]
endpoint = sys.argv[3]

def run_on_endpoint(name, branch):
    import subprocess

    install = subprocess.run(["./install-test.sh "+name+" "+branch], shell=True, encoding="utf_8", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if install.returncode == 1:
        return (install, None)

    result = subprocess.run(["./run-test.sh "+name], shell=True, encoding="utf_8", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return (install, result)

gce = Executor(endpoint_id = endpoint)
future = gce.submit(run_on_endpoint, name, branch)
result = future.result()

os.popen("mkdir -p "+name+"-result").read()
with open(name+"-result/Build.log", "w") as text_file:
    text_file.write("%s" % result[0].stdout)
    text_file.close()
if result[0].returncode == 0:
    with open(name+"-result/Test.log", "w") as text_file:
        text_file.write("%s" % result[1].stdout)
        text_file.close()
