include .env

CR := rindrics
GIT_VERSION := $(shell git describe --tags --always --dirty)
BRANCH := $(shell git symbolic-ref --short HEAD)
VERSION := latest

.PHONY: run
run:
	docker run --rm -it \
	-e PASSWORD=frasyr \
	-v $(PWD):/home/rstudio/frasyr \
	-p 8787:8787 \
	$(CR)/frasyr-dev:$(VERSION)

.PHONY: build
build:
	docker image build \
	-t $(CR)/frasyr-dev:$(BRANCH) \
	-t $(CR)/frasyr-dev:$(GIT_VERSION) \
	-t $(CR)/frasyr-dev:latest \
	.

.PHONY: push
push:
	docker image push $(CR)/frasyr-dev:$(BRANCH)
	docker image push $(CR)/frasyr-dev:$(GIT_VERSION)
	docker image push $(CR)/frasyr-dev:latest
