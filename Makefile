.PHONY: build preview data clean

build:
	@echo "{\"date\": \"$$(gh api /repos/:owner/:repo/commits?per_page=1 --jq '.[0].commit.committer.date' 2>/dev/null || git log -1 --format=%cI)\"}" > data/last_updated.json
	yarn build

preview:
	yarn preview

data:
	@mkdir -p data docs/assets/scenes
	uv run scripts/fetch_scenes.py

clean:
	rm -rf docs/.observable
