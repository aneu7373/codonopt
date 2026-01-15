FROM python:3.11-slim

RUN apt-get update && apt-get install -y viennarna

WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY codonopt/ codonopt/
ENTRYPOINT ["python", "-m", "codonopt.main"]
