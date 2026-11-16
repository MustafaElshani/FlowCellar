"""The Python that assembles pages from the built viewer: reports, the static hub, what serve answers.

Everything here runs without Node and without a browser: a page is text, and these tests read it as text.
"""
from __future__ import annotations

import base64
import gzip
import hashlib
import html
import json
import os
import re
import time
from html.parser import HTMLParser

import pytest
import viewer_synth
from viewer_fake_library import FakeLibrary

from flowcellar.model.schema import schema_errors
from flowcellar.report import api, index_one
from flowcellar.report import html as page
from flowcellar.util.errors import InputError

pytestmark = pytest.mark.viewer
EXEC_TYPES = (None, "", "text/javascript", "module", "application/javascript")


class Page(HTMLParser):
    """The elements of a page that matter for what a browser would run or load."""

    def __init__(self, text: str):
        super().__init__(convert_charrefs=True)
        self.scripts: list[dict] = []
        self.handlers: list[tuple[str, str]] = []
        self.tags: list[str] = []
        self.links: list[dict] = []
        self.meta: dict[str, str] = {}
        self.title = ""
        self.styles: list[str] = []
        self._open: dict | None = None
        self._in: str | None = None
        self.app_text: list[str] = []
        self._depth_app = 0
        self.feed(text)
        self.close()

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        self.tags.append(tag)
        for k, v in attrs:
            if k.lower().startswith("on") or (k.lower() in ("href", "src", "action") and str(v or "").strip().lower().startswith("javascript:")):
                self.handlers.append((tag, k))
        if tag == "script":
            self._open = {"type": a.get("type"), "id": a.get("id"), "src": a.get("src"), "sha": a.get("data-sha256"), "text": ""}
            self._in = "script"
        elif tag == "style":
            self._in = "style"
            self.styles.append("")
        elif tag == "title":
            self._in = "title"
        elif tag == "link":
            self.links.append(a)
        elif tag == "meta" and (a.get("http-equiv") or a.get("name")):
            self.meta[(a.get("http-equiv") or a.get("name")).lower()] = a.get("content") or ""
        if tag == "div" and a.get("id") == "app":
            self._depth_app = 1
        elif self._depth_app and tag == "div":
            self._depth_app += 1

    def handle_endtag(self, tag):
        if tag == "script" and self._open is not None:
            self.scripts.append(self._open)
            self._open = None
        if tag in ("script", "style", "title"):
            self._in = None
        if tag == "div" and self._depth_app:
            self._depth_app -= 1

    def handle_data(self, data):
        if self._in == "script" and self._open is not None:
            self._open["text"] += data
        elif self._in == "style":
            self.styles[-1] += data
        elif self._in == "title":
            self.title += data
        elif self._depth_app:
            self.app_text.append(data)

    @property
    def executable(self):
        return [s for s in self.scripts if s["type"] in EXEC_TYPES]

    @property
    def summary(self) -> str:
        return re.sub(r"\s+", " ", " ".join(self.app_text)).strip()


def sha_csp(text: str) -> str:
    return "'sha256-" + base64.b64encode(hashlib.sha256(text.encode("utf-8")).digest()).decode() + "'"


def bundle_of(made, key):
    return made["bundles"][[s.key for s in made["specs"]].index(key)]


# ---------------------------------------------------------------------------------------------------------
# the built viewer and the template
# ---------------------------------------------------------------------------------------------------------

def test_manifest_matches_the_files_and_the_policy_hashes(viewer_built):
    m = api.viewer_manifest()
    assert m["format"] == "flowcellar.viewer" and m["placeholders"] == list(page.PLACEHOLDERS)
    for name, text in (("viewer.core.js", viewer_built.js), ("viewer.css", viewer_built.css), ("template.html", viewer_built.template)):
        entry = m["files"][name]
        assert entry["sha256"] == hashlib.sha256(text.encode("utf-8")).hexdigest()
        assert entry["bytes"] == len(text.encode("utf-8"))
        assert entry["csp"] == sha_csp(text)
    assert viewer_built.js_csp == m["files"]["viewer.core.js"]["csp"] and viewer_built.css_csp == m["files"]["viewer.css"]["csp"]


def test_template_holds_each_placeholder_exactly_once_in_order(viewer_built):
    found = re.findall(r"%%[A-Z]+%%", viewer_built.template)
    assert found == list(page.PLACEHOLDERS)


def test_fill_is_one_pass_so_a_placeholder_inside_a_value_stays_text(viewer_built):
    text = page.fill(viewer_built, mode="hub", title="%%JS%% %%PAYLOAD%% <b>", summary=page.esc("%%CSS%% & %%BOOT%%"), boot=page.boot_json("hub"))
    p = Page(text)
    assert p.title == "%%JS%% %%PAYLOAD%% <b>"
    assert "%%CSS%% & %%BOOT%%" in p.summary
    assert len(p.executable) == 1 and p.executable[0]["text"] == viewer_built.js
    assert p.styles == [viewer_built.css]


def test_the_built_viewer_holds_nothing_it_must_not(viewer_built):
    js, css = viewer_built.js, viewer_built.css
    assert not re.search(r"</script", js, re.I) and not re.search(r"</style", css, re.I)
    for pat in (r"\bfetch\s*\(", r"XMLHttpRequest", r"new\s+Worker", r"WebSocket", r"sendBeacon", r"\beval\s*\("):
        assert not re.search(pat, js), pat
    for word in ("label missing", "imported, approximate", "sampled"):
        assert word not in js.lower(), word
    urls = [u for u in re.findall(r"https?://[^\s\"'`)]+", js) if not u.startswith("http://www.w3.org/")]
    assert urls == []
    assert not re.search(r"url\(\s*['\"]?(https?:)?//", css) and "@import" not in css


def test_viewer_manifest_says_how_to_build_when_the_viewer_is_missing(monkeypatch, tmp_path):
    from importlib import resources

    class Empty:
        def joinpath(self, *parts):
            return resources.files("flowcellar").joinpath("no_such_folder", *parts)

    monkeypatch.setattr(page, "_dist", lambda: Empty())
    with pytest.raises(InputError) as e:
        api.viewer_manifest()
    assert "viewer is not built" in e.value.message and "build.mjs" in e.value.hint
    with pytest.raises(InputError):
        api.write_report({"bundle_id": "TST00001_0a1b2c3d", "status": {}}, str(tmp_path / "x.html"))
    assert not (tmp_path / "x.html").exists()


# ---------------------------------------------------------------------------------------------------------
# the one-run report
# ---------------------------------------------------------------------------------------------------------

def test_the_synthetic_bundles_are_valid(synth):
    for b in synth["bundles"]:
        assert schema_errors(b) == [], b["bundle_id"]


def test_report_is_one_self_contained_file_with_its_data_inside(viewer_built, synth, tmp_path):
    key = "TST00001_0a1b2c3d"
    bundle = bundle_of(synth, key)
    out = api.write_report(bundle, str(tmp_path / (key + ".flowcellar.html")))
    assert out == str(tmp_path / (key + ".flowcellar.html"))
    raw = open(out, "rb").read()
    assert b"\r\n" not in raw
    text = raw.decode("utf-8")
    p = Page(text)
    # one script runs: the viewer; its hash and the stylesheet's are the policy
    assert len(p.executable) == 1 and p.executable[0]["src"] is None and p.executable[0]["text"] == viewer_built.js
    csp = p.meta["content-security-policy"]
    assert csp == (f"default-src 'none'; script-src {sha_csp(viewer_built.js)}; style-src {sha_csp(viewer_built.css)}; img-src data: blob:; "
                   "font-src 'none'; connect-src 'none'; base-uri 'none'; form-action 'none'")
    assert "'self'" not in csp and "unsafe" not in csp
    # nothing is loaded from anywhere: no src, no stylesheet link, the icon is a data URI
    assert all(s["src"] is None for s in p.scripts)
    assert [link.get("rel") for link in p.links] == ["icon"] and p.links[0]["href"].startswith("data:image/svg+xml;base64,")
    assert not re.search(r"""(?:src|href)\s*=\s*["']?(?:https?:)?//""", text.replace(viewer_built.js, "").replace(viewer_built.css, ""))
    # boot block, index and bundle, in this order, before the script
    by_id = {s["id"]: s for s in p.scripts if s["id"]}
    assert json.loads(by_id["fc-boot"]["text"]) == {"mode": "report", "generator": json.loads(by_id["fc-boot"]["text"])["generator"], "schema": "0.1.1", "run": key}
    assert by_id["fci"]["type"] == "application/vnd.flowcellar.index+gzip"
    data = by_id["fcb:" + key]
    assert data["type"] == "application/vnd.flowcellar.bundle+gzip"
    for s in (by_id["fci"], data):
        assert re.fullmatch(r"[A-Za-z0-9+/]*={0,2}", s["text"])
    gz = base64.b64decode(data["text"])
    assert data["sha"] == hashlib.sha256(gz).hexdigest()
    assert json.loads(gzip.decompress(gz)) == json.loads(json.dumps(bundle))
    index = json.loads(gzip.decompress(base64.b64decode(by_id["fci"]["text"])))
    assert index["format"] == "flowcellar.index" and index["runs"]["n"] == 1
    order = [text.index('id="app"'), text.index('id="fc-boot"'), text.index('id="fci"'), text.index('id="fcb:'), text.index("<script>")]
    assert order == sorted(order)


def test_report_text_summary_reads_without_scripts(viewer_built, synth, tmp_path):
    bundle = bundle_of(synth, "TST00002_1b2c3d4e")
    p = Page(open(api.write_report(bundle, str(tmp_path / "r.html")), encoding="utf-8").read())
    s = p.summary
    st = bundle["status"]
    assert "FlowCellar run report" in s and "TST00002" in s
    assert "Warning" in s and st["headline"] in s
    assert "Run health" in s and "Read data" in s
    for rule in st["rules"]:
        if rule["outcome"] == "triggered":
            assert rule["message"] in s
    assert "Pores at the first scan" in s and f"{bundle['modules']['fc']['summary']['fc_pores_first_scan']:,}" in s
    assert "Reads" in s and "Estimated bases" in s
    assert p.title == "TST00002 · Warning · FlowCellar run report"


def test_report_of_a_run_in_quick_detail_says_so_far(viewer_built, synth, tmp_path):
    p = Page(open(api.write_report(bundle_of(synth, "TST00004_3d4e5f6a"), str(tmp_path / "r.html")), encoding="utf-8").read())
    assert p.title == "TST00004 · Pass so far · FlowCellar run report" and "Pass so far" in p.summary


def test_report_of_a_folder_with_raw_data_only(viewer_built, synth, tmp_path):
    p = Page(open(api.write_report(bundle_of(synth, "TST00006_7b8c9d0e"), str(tmp_path / "r.html")), encoding="utf-8").read())
    assert "Incomplete" in p.summary and "Only raw data was found" in p.summary


def test_names_that_try_to_be_markup_stay_text_in_a_report(viewer_built, synth, tmp_path):
    key = "TST00009_0e1f2a3b"
    bundle = bundle_of(synth, key)
    assert bundle["identity"]["protocol_run"]["sample_id"] == viewer_synth.HOSTILE_RUN_NAME
    text = open(api.write_report(bundle, str(tmp_path / "r.html")), encoding="utf-8").read()
    p = Page(text)
    # exactly the four script elements a report has, and only one of them runs
    assert sorted(s["id"] or "viewer" for s in p.scripts) == ["fc-boot", "fcb:" + key, "fci", "viewer"]
    assert len(p.executable) == 1 and p.executable[0]["text"] == viewer_built.js
    assert p.handlers == []
    assert "img" not in p.tags and "svg" not in p.tags and "iframe" not in p.tags and "object" not in p.tags
    # the name is there, as text
    assert viewer_synth.HOSTILE_RUN_NAME in p.summary
    outside = text.replace(viewer_built.js, "").replace(viewer_built.css, "")
    assert "window.__pwned" not in re.sub(r"&lt;.*?&gt;", "", outside.replace(html.escape(viewer_synth.HOSTILE_RUN_NAME), ""))
    assert "<script>window.__pwned" not in outside and "onerror=" not in outside.replace(html.escape(viewer_synth.HOSTILE_RUN_NAME), "")


@pytest.mark.parametrize("name", ["</script><script>alert(1)</script>", "--><script>alert(1)</script><!--", "]]></script>", "%%JS%%", " x ", "' onload='x", '"><svg/onload=1>'])
def test_any_run_name_leaves_the_page_structure_alone(viewer_built, synth, tmp_path, name):
    bundle = json.loads(json.dumps(bundle_of(synth, "TST00003_2c3d4e5f")))
    bundle["identity"]["protocol_run"]["sample_id"] = name
    bundle["status"]["headline"] = name
    bundle["identity"]["flow_cell"]["id"] = "TST00003"
    p = Page(open(api.write_report(bundle, str(tmp_path / "r.html")), encoding="utf-8").read())
    assert len(p.scripts) == 4 and len(p.executable) == 1 and p.executable[0]["text"] == viewer_built.js
    assert p.handlers == [] and "svg" not in p.tags
    assert name.strip() in p.summary or re.sub(r"\s+", " ", name).strip() in p.summary


def test_boot_block_never_holds_a_markup_character():
    assert "<" not in page.boot_json("report", "TST00001_0a1b2c3d") and ">" not in page.boot_json("hub")
    with pytest.raises(InputError):
        page.boot_json("report", "</script><script>")
    with pytest.raises(InputError):
        page.run_sidecar('x");alert(1);//', b"")
    with pytest.raises(InputError):
        page.data_element(page.BUNDLE_TYPE, 'fcb:x" onload="y', b"")


def test_report_refuses_a_bundle_without_a_verdict_or_a_run_key(viewer_built, synth, tmp_path):
    b = dict(bundle_of(synth, "TST00001_0a1b2c3d"))
    del b["status"]
    with pytest.raises(InputError):
        api.write_report(b, str(tmp_path / "r.html"))
    b2 = dict(bundle_of(synth, "TST00001_0a1b2c3d"), bundle_id="../../etc/passwd")
    with pytest.raises(InputError):
        api.write_report(b2, str(tmp_path / "r.html"))
    assert not (tmp_path / "r.html").exists()


def test_report_is_written_whole_and_twice_the_same(viewer_built, synth, tmp_path):
    b = bundle_of(synth, "TST00003_2c3d4e5f")
    index = index_one.plain_index(b)
    one = open(api.write_report(b, str(tmp_path / "a.html"), index=index), "rb").read()
    two = open(api.write_report(b, str(tmp_path / "b.html"), index=index), "rb").read()
    assert one == two
    assert [f for f in os.listdir(tmp_path) if f.endswith(".tmp")] == []


def test_plain_index_has_the_shape_of_the_library_index(synth):
    b = bundle_of(synth, "TST00001_0a1b2c3d")
    ix = index_one.plain_index(b)
    for member in ("format", "schema", "written", "meta", "library", "runs", "samples", "projects", "roots", "activity", "baselines", "facets", "notices", "progress"):
        assert member in ix
    cols = ix["runs"]["columns"]
    assert cols["id"] == ["TST00001_0a1b2c3d"] and cols["verdict"] == [b["status"]["verdict"]] and cols["headline"] == [b["status"]["headline"]]
    assert cols["health"] == ["pass"] and cols["read_data"] == ["complete"] and cols["reads"] == [b["modules"]["yield"]["summary"]["yield_reads"]]
    assert len(ix["samples"]) == 4


def test_slice_index_keeps_one_row_and_nothing_of_the_other_runs(synth):
    ix = index_one.slice_index(synth["index"], "TST00002_1b2c3d4e")
    assert ix["runs"]["n"] == 1 and ix["runs"]["columns"]["id"] == ["TST00002_1b2c3d4e"]
    assert {s["run"] for s in ix["samples"]} == {"TST00002_1b2c3d4e"}
    assert ix["roots"] == [] and ix["library"]["path"] is None
    text = json.dumps(ix)
    for other in ("TST00001", "TST00009", "/data/sequencer", "Demo Beta"):
        assert other not in text
    assert index_one.slice_index(synth["index"], "TST09999_00000000") is None


# ---------------------------------------------------------------------------------------------------------
# the static hub
# ---------------------------------------------------------------------------------------------------------

def test_hub_page_holds_the_viewer_and_no_data(viewer_built, library):
    lib, made = library
    res = api.write_hub(lib)
    assert res["index_html"] == os.path.join(lib.hub_dir(), "index.html")
    text = open(res["index_html"], encoding="utf-8").read()
    p = Page(text)
    assert [s["id"] for s in p.scripts if s["id"]] == ["fc-boot"]
    assert json.loads(p.scripts[0]["text"])["mode"] == "hub" and json.loads(p.scripts[0]["text"])["run"] is None
    assert len(p.executable) == 1 and p.executable[0]["text"] == viewer_built.js
    assert p.meta["content-security-policy"] == (
        f"default-src 'none'; script-src {sha_csp(viewer_built.js)} 'self'; style-src {sha_csp(viewer_built.css)}; img-src data: blob:; "
        "font-src 'none'; connect-src 'none'; base-uri 'none'; form-action 'none'")
    # no name of the library is in the page: it changes only when the viewer changes
    for word in ("TST0000", "Demo", "/data/"):
        assert word not in text.replace(viewer_built.js, "").replace(viewer_built.css, "")
    assert "needs a browser with scripts" in p.summary


def test_hub_sidecars_are_exactly_one_call_each(viewer_built, library):
    lib, made = library
    api.write_hub(lib)
    hub = lib.hub_dir()
    text = open(os.path.join(hub, "index.js"), encoding="ascii").read()
    m = re.fullmatch(r'__FC_INDEX__\("([A-Za-z0-9+/]*={0,2})"\);\n', text)
    assert m and json.loads(gzip.decompress(base64.b64decode(m.group(1)))) == lib.index()
    assert base64.b64decode(m.group(1)) == lib.index_gzip()
    assert sorted(os.listdir(os.path.join(hub, "data"))) == [k + ".js" for k in lib.run_keys()]
    for key in lib.run_keys():
        t = open(os.path.join(hub, "data", key + ".js"), encoding="ascii").read()
        m = re.fullmatch(r'__FCB_SIDECAR__\("([A-Z0-9]+_[0-9a-f]+)","([A-Za-z0-9+/]*={0,2})","([0-9a-f]{64})"\);\n', t)
        assert m and m.group(1) == key
        gz = base64.b64decode(m.group(2))
        assert hashlib.sha256(gz).hexdigest() == m.group(3)
        b = json.loads(gzip.decompress(gz))
        assert b["bundle_id"] == key and "status" in b and "relations" in b and "user" in b


def _mtimes(hub):
    out = {"index.html": os.stat(os.path.join(hub, "index.html")).st_mtime_ns, "index.js": os.stat(os.path.join(hub, "index.js")).st_mtime_ns}
    for f in os.listdir(os.path.join(hub, "data")):
        out[f] = os.stat(os.path.join(hub, "data", f)).st_mtime_ns
    return out


def test_hub_rewrites_only_the_sidecars_of_changed_runs(viewer_built, library):
    lib, made = library
    hub = lib.hub_dir()
    first = api.write_hub(lib)
    n = len(lib.run_keys())
    assert len(first["written"]) == n + 2 and first["removed"] == []
    before = _mtimes(hub)
    time.sleep(0.02)
    second = api.write_hub(lib)
    assert second["written"] == [os.path.join(hub, "index.js")]          # the index always, nothing else
    after = _mtimes(hub)
    assert after["index.html"] == before["index.html"]
    assert all(after[f] == before[f] for f in before if f.endswith(".js") and f != "index.js")

    # a run is judged again: its status file is newer, and only its sidecar is written
    key = "TST00002_1b2c3d4e"
    sf_path = os.path.join(lib.path, "status", key + ".status.json")
    sf = json.load(open(sf_path))
    sf["status"]["headline"] = "A new sentence."
    time.sleep(0.02)
    with open(sf_path, "w") as fh:
        json.dump(sf, fh)
    third = api.write_hub(lib)
    assert sorted(third["written"]) == sorted([os.path.join(hub, "index.js"), os.path.join(hub, "data", key + ".js")])
    t = open(os.path.join(hub, "data", key + ".js")).read()
    assert json.loads(gzip.decompress(base64.b64decode(t.split('"')[3])))["status"]["headline"] == "A new sentence."

    # the caller names a run: it is written whether or not it changed
    fourth = api.write_hub(lib, only=["TST00003_2c3d4e5f", "TST09999_00000000"])
    assert sorted(fourth["written"]) == sorted([os.path.join(hub, "index.js"), os.path.join(hub, "data", "TST00003_2c3d4e5f.js")])

    # a run leaves the library: its sidecar goes
    os.unlink(os.path.join(lib.path, "bundles", "TST00007_8c9d0e1f.fcb.json.gz"))
    os.unlink(os.path.join(lib.path, "status", "TST00007_8c9d0e1f.status.json"))
    fifth = api.write_hub(lib)
    assert fifth["removed"] == [os.path.join(hub, "data", "TST00007_8c9d0e1f.js")]
    assert not os.path.exists(fifth["removed"][0])
    assert [f for f in os.listdir(os.path.join(hub, "data")) if f.endswith(".tmp")] == []


def test_hub_without_sidecars_leaves_no_stale_data_behind(viewer_built, library):
    lib, made = library
    hub = lib.hub_dir()
    api.write_hub(lib)
    key = "TST00001_0a1b2c3d"
    time.sleep(0.02)
    os.utime(os.path.join(lib.path, "status", key + ".status.json"))
    res = api.write_hub(lib, sidecars=False)
    assert res["removed"] == [os.path.join(hub, "data", key + ".js")]
    assert res["written"] == [os.path.join(hub, "index.js")]
    assert len(os.listdir(os.path.join(hub, "data"))) == len(lib.run_keys()) - 1


def test_hub_can_write_one_report_per_run(viewer_built, library):
    lib, made = library
    res = api.write_hub(lib, reports=True)
    reports = sorted(os.listdir(os.path.join(lib.path, "reports")))
    assert reports == [k + ".flowcellar.html" for k in lib.run_keys()]
    assert sum(1 for p in res["written"] if p.endswith(".flowcellar.html")) == len(reports)
    p = Page(open(os.path.join(lib.path, "reports", "TST00001_0a1b2c3d.flowcellar.html"), encoding="utf-8").read())
    index = json.loads(gzip.decompress(base64.b64decode(next(s for s in p.scripts if s["id"] == "fci")["text"])))
    assert index["runs"]["n"] == 1 and index["runs"]["columns"]["id"] == ["TST00001_0a1b2c3d"]
    again = api.write_hub(lib, reports=True)
    assert not any(p.endswith(".flowcellar.html") for p in again["written"])
    # a report made with another build of the viewer is written again
    old = os.path.join(lib.path, "reports", "TST00003_2c3d4e5f.flowcellar.html")
    text = open(old, encoding="utf-8").read()
    with open(old, "w", encoding="utf-8", newline="") as fh:
        fh.write(re.sub(r"script-src &#x27;sha256-[A-Za-z0-9+/=]+&#x27;", "script-src &#x27;sha256-" + "B" * 43 + "=&#x27;", text, count=1))
    third = api.write_hub(lib, reports=True)
    assert [p for p in third["written"] if p.endswith(".flowcellar.html")] == [old]
    assert open(old, encoding="utf-8").read() == text


def test_hub_command_prints_where_the_hub_is(viewer_built, library, monkeypatch):
    from click.testing import CliRunner

    from flowcellar.cli import cli
    from flowcellar.library import api as library_api

    lib, made = library
    monkeypatch.setattr(library_api, "open_library", lambda path, **kw: FakeLibrary(path))
    r = CliRunner().invoke(cli, ["hub", lib.path])
    assert r.exit_code == 0, r.output
    assert r.output.startswith("hub ready: " + os.path.join(lib.path, "hub", "index.html"))
    assert f"{len(lib.run_keys())} runs" in r.output and f"{len(lib.run_keys())} data files written" in r.output
    r = CliRunner().invoke(cli, ["hub", lib.path, "--no-sidecars"])
    assert r.exit_code == 0 and "open it with 'flowcellar serve'" in r.output


# ---------------------------------------------------------------------------------------------------------
# what serve answers
# ---------------------------------------------------------------------------------------------------------

def test_serve_answers_four_addresses_and_nothing_else(viewer_built, library):
    lib, made = library
    key = "TST00001_0a1b2c3d"
    # no hub was written yet: there is no page to serve (the serve command then says how to write one) ...
    assert api.serve_file(lib, "/") is None and api.serve_file(lib, "/index.html") is None
    # ... but the index and a run's data can be made from the library, the same bytes a hub would hold
    side = api.serve_file(lib, f"/data/{key}.js")
    assert not os.path.exists(os.path.join(lib.hub_dir(), "data"))
    api.write_hub(lib)
    hub = lib.hub_dir()
    root = api.serve_file(lib, "/")
    assert root.status == 200 and root.content_type == "text/html; charset=utf-8"
    assert root.body == open(os.path.join(hub, "index.html"), "rb").read() == api.serve_file(lib, "/index.html").body
    p = Page(root.body.decode("utf-8"))
    for answer in (root, api.serve_file(lib, "/index.js"), api.serve_file(lib, f"/data/{key}.js")):
        assert answer.headers["Content-Security-Policy"] == p.meta["content-security-policy"]
        assert answer.headers["X-Content-Type-Options"] == "nosniff" and answer.headers["Referrer-Policy"] == "no-referrer"
        assert answer.headers["Cross-Origin-Opener-Policy"] == "same-origin" and answer.headers["Cache-Control"] == "no-store"
    ix = api.serve_file(lib, "/index.js?t=3")
    assert ix.content_type == "text/javascript; charset=utf-8" and ix.body == open(os.path.join(hub, "index.js"), "rb").read()
    assert side.body == open(os.path.join(hub, "data", key + ".js"), "rb").read() == api.serve_file(lib, f"/data/{key}.js").body
    for path in ("/data/", "/data/TST09999_00000000.js", "/data/../LIBRARY", "/data/..%2f..%2fLIBRARY", "/LIBRARY", "/flowcellar.toml", "/bundles/" + key + ".fcb.json.gz",
                 "/favicon.ico", "/data/" + key, "/data/" + key + ".js/", "//index.js", "/index.js/../LIBRARY", "/data/" + key + ".JS", "", "/hub/index.html",
                 "/index.html/", "/data//index.html", "/index.next.js", "/reports/" + key + ".flowcellar.html"):
        assert api.serve_file(lib, path) is None, path


def test_serve_states_the_policy_of_the_page_it_serves(viewer_built, library):
    """A hub page written by an older viewer carries that viewer's script hash; the header must say the same."""
    lib, made = library
    api.write_hub(lib)
    index_html = os.path.join(lib.hub_dir(), "index.html")
    text = open(index_html, encoding="utf-8").read()
    old = "sha256-" + "A" * 43 + "="
    stale = re.sub(r"script-src &#x27;sha256-[A-Za-z0-9+/=]+&#x27;", "script-src &#x27;" + old + "&#x27;", text, count=1)
    assert stale != text
    with open(index_html, "w", encoding="utf-8", newline="") as fh:
        fh.write(stale)
    served = api.serve_file(lib, "/")
    assert served.body.decode("utf-8") == stale
    assert f"script-src '{old}' 'self'" in served.headers["Content-Security-Policy"]
    assert served.headers["Content-Security-Policy"] == Page(stale).meta["content-security-policy"]


def test_serve_never_hands_out_a_stale_sidecar(viewer_built, library):
    lib, made = library
    api.write_hub(lib)
    key = "TST00002_1b2c3d4e"
    sf_path = os.path.join(lib.path, "status", key + ".status.json")
    sf = json.load(open(sf_path))
    sf["status"]["headline"] = "Judged again."
    time.sleep(0.02)
    with open(sf_path, "w") as fh:
        json.dump(sf, fh)
    body = api.serve_file(lib, f"/data/{key}.js").body.decode()
    assert json.loads(gzip.decompress(base64.b64decode(body.split('"')[3])))["status"]["headline"] == "Judged again."
