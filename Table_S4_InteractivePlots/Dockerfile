# Use an official R base image
FROM r-base:4.3.1

# Set the working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install required R packages and their dependencies
RUN R -e "install.packages(c('shiny', 'plotly', 'dplyr', 'DT'))"

# Copy required files into the container
COPY . /app/

# Run the R Shiny app
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
