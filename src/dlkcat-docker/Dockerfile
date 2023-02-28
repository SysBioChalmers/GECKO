# syntax=docker/dockerfile:1

FROM python:3.9-slim

WORKDIR ./

COPY requirements.txt ./

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD [ "python", "./DLKcat.py" ]


