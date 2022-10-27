.PHONY: conf
conf:
	cmake -S . -B build

.PHONY: build
build:
	cmake --build build

.PHONY: test
test:
	cmake --build build --target test

.PHONY: run
run:
	./$(pwd)/build/apps/app

.PHONY: brun
brun:
	cmake -S . -B build >> build/buildlog.txt
	./$(pwd)/build/apps/app
