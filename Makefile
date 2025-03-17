# Let's dynamically generate our list of containers to build.
# Find all directories that begin with lowercase characters
# exclude .git and .
# then read from that generated file.
# and build/push em to dockerhub?

CONTAINERS_FILE = containers.txt

CONTAINERS = $(shell cat $(CONTAINERS_FILE) 2>/dev/null | grep -v "^$$" | tr -d ' \t\r')

generate-containers-list:
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

clean-containers:
	rm -f $(CONTAINERS_FILE)

# force regen the container list
refresh-modules: clean-containers generate-containers-list 
