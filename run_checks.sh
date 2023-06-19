# Script to run GitHub Actions checks locally.
# This isn't a 1:1 copy of the GitHub workflow (e.g. only runs on one OS & Python 
# version) but it's generally close enough to catch most issues before pushing.

echo -e "RUNNING CHECKS (`date`)" && \
echo -e "\nRunning black..." && black . --check && \
echo -e "\nRunning isort..." && isort . --check && \
echo -e "\nRunning pylint..." && pylint $(git ls-files '*.py') && \
echo -e "\nALL DONE"
