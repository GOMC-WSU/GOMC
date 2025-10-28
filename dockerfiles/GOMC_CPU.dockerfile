# Use the Debian base image
FROM debian:bullseye

# Install necessary dependencies
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    git \
    wget \
    python \
    cmake

# Clone the GOMC repository
RUN git clone -b wolf_fitz --single-branch https://github.com/GregorySchwing/GOMC.git /gomc

# Set the working directory to the GOMC repository
WORKDIR /gomc

# Build the code using metamake.sh
RUN ./metamake.sh

# Add /gomc/bin to the PATH
ENV PATH="/gomc/bin:${PATH}"

# Command to run when the container starts
CMD ["bash"]
