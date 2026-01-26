# Use a slim Python 3.11 base
FROM python:3.11-slim

WORKDIR /app

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the whole codonopt package
COPY codonopt/ codonopt/

# Copy optional README
COPY README.md .

# Entrypoint points to main.py
ENTRYPOINT ["python", "-m", "codonopt.main"]