## TODO: fix this doc and add more functionality for template handling, linting, validation, etc.

# add simple argparsing to the makefile
CMD := $(word 1,$(MAKECMDGOALS))
ARG := $(filter-out $(CMD),$(MAKECMDGOALS))

# don't error on missing targets?
%::
	@true

.PHONY: new-deb
new-deb:
	@if [ -z "$(ARG)" ]; then echo "ERROR: container name required. Try again."; exit 1; fi
	@NEWIMG="$(ARG)" ; \
	TMPS="templates" ; \
	MAKEFILE="templates/Makefile-template" ; \
	DOCKERFILE="templates/Dockerfile-bookworm-slim" ; \
	echo "----- Initializing new debian bookworm container: $$NEWIMG -----" \
	&& mkdir -p "$$NEWIMG" \
	&& sed "s|IMGNAMEPLACEHOLDER|$$NEWIMG|g" "$$DOCKERFILE" >"$$NEWIMG/Dockerfile" \
	&& sed "s|IMGNAMEPLACEHOLDER|$$NEWIMG|g" "$$TMPS/Makefile-template" >"$$NEWIMG/Makefile"

.PHONY: new-py
new-py:
	@if [ -z "$(ARG)" ]; then echo "ERROR: container name required. Try again."; exit 1; fi
	@NEWIMG="$(ARG)" ; \
	MAKEFILE="templates/Makefile-template"; \
	DOCKERFILE="templates/Dockerfile-uv-python3-13" ; \
	echo "----- Initializing new python3.13 container: $$NEWIMG -----" \
	&& mkdir -p "$$NEWIMG" \
	&& sed "s|IMGNAMEPLACEHOLDER|$$NEWIMG|g" "$$DOCKERFILE" >"$$NEWIMG/Dockerfile" \
	&& sed "s|IMGNAMEPLACEHOLDER|$$NEWIMG|g" "$$MAKEFILE" >"$$NEWIMG/Makefile"

.PHONY: new-cuda
new-cuda:
	@if [ -z "$(ARG)" ]; then echo "ERROR: container name required. Try again."; exit 1; fi
	@NEWIMG="$(ARG)" ; \
	TMPS="templates" ; \
	MAKEFILE="templates/Makefile-template" ; \
	DOCKERFILE="templates/Dockerfile-cuda13-ubuntu24-04" ; \
	echo "----- Initializing new debian bookworm container: $$NEWIMG -----" \
	&& mkdir -p "$$NEWIMG" \
	&& sed "s|IMGNAMEPLACEHOLDER|$$NEWIMG|g" "$$DOCKERFILE" >"$$NEWIMG/Dockerfile" \
	&& sed "s|IMGNAMEPLACEHOLDER|$$NEWIMG|g" "$$TMPS/Makefile-template" >"$$NEWIMG/Makefile"

clean-containers:
	rm $(CONTAINERS_FILE)

# force regen the container list
refresh-modules: clean-containers generate-containers-list

# Let's dynamically generate our list of containers to build.
generate-containers-list:
	CONTAINERS_FILE = containers.txt
	CONTAINERS = $(shell cat $(CONTAINERS_FILE) 2>/dev/null | grep -v "^$$" | tr -d ' \t\r')
	find . -maxdepth 1 -type d -name "[a-z]*" -not -name ".git" -printf "%f\n" |\
		sort |\
		uniq > $(CONTAINERS_FILE)
	@echo "Generated list of containers to build in $(CONTAINERS_FILE)"

all:
	for dir in $(CONTAINERS); do \
		echo begin $$dir; \
		(cd $$dir; ${MAKE}); \
		echo end $$dir; \
		echo ==========; \
	done
