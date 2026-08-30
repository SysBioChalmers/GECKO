Use this command to obtain the multiarch image:

```
docker buildx build --platform linux/amd64,linux/arm64/v8 --output type=image,name=ghcr.io/sysbiochalmers/dlkcat-gecko,push=true  -t ghcr.io/sysbiochalmers/dlkcat-gecko:0.1-multiarch . 
```
