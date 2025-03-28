#!/bin/bash

# Build image
docker build -t trplet-app .

# Run container
docker run --rm -ti -p 3838:3838 trplet-app

