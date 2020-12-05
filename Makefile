.PHONY: all lint test test-cov viz-summarize install dev clean distclean

PYTHON ?= python

all: viz-summarize

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	py.test --cov=q2_demux

q2_demux/_summarize/assets/dist:
	cd q2_demux/_summarize/assets && \
	npm install --no-save && \
	npm run build && \
	cp licenses/* dist/

viz-summarize: q2_demux/_summarize/assets/dist

install: all
	$(PYTHON) setup.py install

dev: all
	pip install -e .

clean: distclean
	rm -rf q2_demux/_summarize/assets/node_modules

distclean:
	rm -rf q2_demux/_summarize/assets/dist
