docker build --platform linux/amd64 -t cando .
docker run -p 8888:8888 --platform linux/amd64 -it cando

