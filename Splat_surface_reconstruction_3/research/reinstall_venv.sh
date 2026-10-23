# Source - https://stackoverflow.com/a/61080118
# Posted by Steve Rossiter
# Retrieved 2026-10-11, License - CC BY-SA 4.0

uv pip freeze > to-uninstall.txt
uv pip uninstall -r to-uninstall.txt
uv pip install -r requirements.txt
# delete the temporary file containing the list of packages to uninstall
rm to-uninstall.txt