# Packaging & Release — Alpha distribution plan

**Status:** Foundation IMPLEMENTED (2026-06-18). Local macOS `.dmg` built + validated; CI authored.
**Date:** 2026-06-16 (design) · 2026-06-18 (implementation)
**Scope:** How to ship the Rust `navigator` binary as installable images for macOS
(Apple Silicon + Intel), Linux (x86-64 + ARM), and Windows (x86-64 + ARM), and how to
handle the bundled ancestry assets.

---

## Implementation status (2026-06-18)

**Landed (committed on rust-rewrite):**
- Deleted the legacy Scala `jpackage` `release.yml`.
- `.cargo/config.toml` — `target-cpu=ivybridge` for the three x86-64 triples (ARM untouched).
- `[package.metadata.packager]` in `crates/navigator-ui/Cargo.toml` (product `DUNavigator`,
  id `com.decodingus.navigator`, macOS min 11.0, icons, `before-packaging-command` → staging,
  `resources` → bundled `ancestry/`). Placeholder icon at `crates/navigator-ui/icons/`.
- `packaging/stage-assets.sh` stages the full Option-A bundle (copies from `~/.decodingus/ancestry`
  locally; CDN path is a TODO). `packaging/staging/` is git-ignored.
- First-run **asset seeding**: `navigator_app::seed_bundled_assets` / `seed_assets_from` (copies
  missing bundled assets → `~/.decodingus/ancestry/`, never overwriting a CDN-refreshed file),
  called at `main()` startup (both GUI + headless). Unit-tested.
- New cargo-packager `release.yml` (matrix: macOS aarch64/x86_64, Linux x86_64/aarch64, Windows
  x86_64; collect + `SHA256SUMS` + gh-release, prerelease for `-alpha`).

**Validated locally (macOS arm64):** `cargo packager --release -p navigator-ui -f dmg` →
`DUNavigator_0.1.0_aarch64.dmg` with `Contents/Resources/ancestry/` fully populated (incl. the
102 MB IBD panel), `icon.icns`, `CFBundleIdentifier=com.decodingus.navigator`,
`LSMinimumSystemVersion=11.0`; the bundled binary runs; seeding into a fresh
`NAVIGATOR_REFGENOME_DIR` confirmed (11 assets).

**macOS universal2 (2026-06-24, DONE in CI + validated locally):** the CI now ships **one
universal `.dmg`** instead of per-arch ones. cargo-packager only *packages* a binary (it does not
build a universal one), so the `package-macos` job builds both slices, `lipo`s them into
`target/universal-apple-darwin/release/navigator`, then `cargo packager --target
universal-apple-darwin -f dmg`. Validated locally on macos-14: `DUNavigator_0.1.0_universal.dmg`,
the bundled `Contents/MacOS/navigator` is a fat `x86_64 arm64` binary with `minos 11.0` and
`CFBundleIdentifier=com.decodingus.navigator`. The Intel slice gets `target-cpu=ivybridge` from
`.cargo/config.toml`.

**Linux glibc floor (2026-06-24, decided = 2.28 via container; CI authored, NOT yet run):** the
`package-linux` job builds **and packages inside `quay.io/pypa/manylinux_2_28_*`** (AlmaLinux 8,
glibc 2.28) for both x86_64 and aarch64, so the binary *and* the AppImage-bundled GTK/D-Bus libs all
target 2.28 (Ubuntu 20.04+/Debian 10+/RHEL 8+). `cargo-zigbuild` was evaluated and rejected: it only
pins *our* binary's glibc, not the libs linuxdeploy bundles into the AppImage, so it gives no real
floor benefit for a GUI app (a local macOS cross-build also confirmed the Rust + vendored-C deps —
sqlite/bzip2/lzma/ring/noodles — cross-compile cleanly; only the GUI system libs need the container).
manylinux_2_28 is chosen because it runs GitHub's node20 actions (2.28 is node20's floor); AppImage
FUSE is bypassed with `APPIMAGE_EXTRACT_AND_RUN=1`. **Unverified — no local Linux/Docker here; needs
a CI run to settle the exact `dnf` package set + AppImage tooling.**

**General CI fix (2026-06-24):** cargo-packager 0.11.8 only *packages* a pre-built binary — it does
**not** build with `--target` (verified locally: it errors "No such file" instead of compiling). So
every job now runs an explicit `cargo build`/`lipo` before `cargo packager`; the original matrix
(which called `cargo packager --target …` with no prior build) would have failed.

**Windows (2026-06-24, validated via cross-build):** the `package-windows` job builds the
`x86_64-pc-windows-msvc` binary then `cargo packager … --formats nsis` → `*_x64-setup.exe`. Validated
locally from macOS with `cargo-xwin` (clang-cl + Windows SDK + nasm for `ring`): the full
workspace + all vendored-C deps (sqlite/bzip2/lzma/ring) and the egui/eframe/rfd/keyring GUI stack
compile clean, and `cargo packager` produced a working NSIS self-extracting installer (PE32, ~40 MB,
the ~190 MB binary+assets compressed) using a host `makensis`. The CI uses native `windows-latest`
(MSVC) — even closer to target — so this is a strong proxy. ivybridge floor comes from
`.cargo/config.toml`; the collect step ships only `*-setup.exe`, not the bare `navigator.exe`. Alpha
is **unsigned** (cargo-packager warns + skips signing off-Windows); Authenticode / Azure Trusted
Signing comes before beta.

**Still CI-time / unverified here (need a tagged run + secrets):** macOS notarization (`APPLE_*`),
the Linux container job (above — needs a CI run), Windows signing (unsigned for Alpha), and the CI
asset-staging CDN source (`NAVIGATOR_ASSET_SRC`/CDN).

---

## TL;DR / recommendation

1. **Delete the legacy `release.yml`** — it builds the Scala fat JAR with `jpackage`
   and is unrelated to the Rust binary.
2. **Build natively on each OS runner** (no cross-compilation toolchain juggling). The
   stack makes this clean: TLS is rustls (no OpenSSL), and the only C deps (SQLite,
   bzip2, xz, `ring`) are vendored into the binary. Each target is a single
   self-contained executable.
3. **Use `cargo-packager`** (CrabNebula) to produce the OS-native images: `.app`+`.dmg`,
   `.msi`/NSIS, AppImage/`.deb`/`.rpm`. It is GUI-first, covers all three OSes including
   Windows, and supports macOS notarization + Windows Authenticode + an updater.
   - Alternative: `cargo-dist` (still maintained, v0.31, 2026-02) if we want the
     turnkey "generate the whole CI + curl|sh installer + auto-update" story. It is more
     CLI-distribution-oriented but now does `.dmg`/`.msi` and signing too.
4. **CPU floor = Ivy Bridge** is a hard constraint already documented in
   `simd-optimization-targets.md`. Compile the x86-64 slices with
   `-C target-cpu=ivybridge` (= `x86-64-v2` + AVX). **Never** `x86-64-v3` (AVX2/FMA →
   `SIGILL` on the baseline Mac Pro 2013). ARM64 always has NEON; nothing to set.
5. **Bundle all 160 MB of ancestry assets** in the image (full offline installer,
   ≈180–200 MB per platform — *decided*) and seed them to `~/.decodingus/ancestry/` on
   first run, verified against the existing `AssetManifest`.

**Decisions locked (2026-06-16):** full offline installer · macOS notarized only (Windows
unsigned w/ documented work-around for Alpha) · `cargo-packager`. Implementation deferred —
this stays research for now.

---

## 1. Target matrix

| OS | Arch | Rust triple | Image | CPU floor handling |
|----|------|-------------|-------|--------------------|
| macOS | Apple Silicon | `aarch64-apple-darwin` | `.app` in `.dmg` (universal) | NEON baseline |
| macOS | Intel | `x86_64-apple-darwin` | (same universal binary) | `target-cpu=ivybridge` |
| Linux | x86-64 | `x86_64-unknown-linux-gnu` | AppImage + `.deb`/`.rpm` + `.tar.gz` | `target-cpu=ivybridge` |
| Linux | ARM | `aarch64-unknown-linux-gnu` | AppImage + `.deb`/`.rpm` + `.tar.gz` | NEON baseline |
| Windows | x86-64 | `x86_64-pc-windows-msvc` | `.msi` (or NSIS) + portable `.zip` | `target-cpu=ivybridge` |
| Windows | ARM | `aarch64-pc-windows-msvc` | `.msi` + `.zip` (defer for Alpha) | NEON baseline |

**macOS:** ship a **universal2 binary** (lipo of both slices) inside one `.dmg`, so a
single download works on both Macs. Apple `clang` cross-compiles both arches from one
runner; the vendored C deps build for both. Set `MACOSX_DEPLOYMENT_TARGET=11.0` (Big Sur;
the Mac Pro 2013 baseline tops out at Monterey, so this is safe and `aarch64` requires
≥11.0 anyway).

**Linux glibc floor:** GitHub `ubuntu-22.04` = glibc 2.35, which is too new for older
distros. For broad "Ivy Bridge or newer" reach, either (a) build inside an old-glibc
container (e.g. manylinux / Debian 10 → glibc 2.28), or (b) use **`cargo-zigbuild`** with
a pinned glibc (`x86_64-unknown-linux-gnu.2.28`). AppImage does **not** bundle glibc, so
the build-glibc choice still gates compatibility. ARM Linux builds natively on the
`ubuntu-24.04-arm` runner (avoids cross C-toolchain pain entirely).

**Windows ARM:** `aarch64-pc-windows-msvc` cross-builds from an x64 runner but is the
least-exercised slice; recommend **deferring it past Alpha** unless an ARM Windows tester
exists.

---

## 2. Why native-per-runner instead of cross-compiling

The usual reason to cross-compile is to avoid N runners; the usual reason *not* to is
that C/asm deps need a cross C-toolchain and sysroot. Here:

- `libsqlite3-sys` (sqlx), `bzip2-sys` + `lzma-sys` (noodles CRAM), `ring` — all compile C
  or asm via `cc`. Cross-compiling these needs a matching cross toolchain.
- `flate2` uses the pure-Rust `miniz_oxide` backend (no system zlib). ✅
- `reqwest` is `default-features = false, features = ["rustls-tls"]` everywhere — **no
  OpenSSL**. ✅

Native-per-runner means the host compiler already has the right toolchain, so the C deps
"just build." The only places we deliberately cross within a runner are: macOS x86_64
slice (Apple clang handles it) and, optionally, old-glibc Linux via zig. This is the
lowest-friction path and what `cargo-packager`/`cargo-dist` assume by default.

---

## 3. CPU baseline (ties into `simd-optimization-targets.md`)

The reference machine is a **Mac Pro 2013 (Xeon E5 v2 = Ivy Bridge-EP)**. Ivy Bridge has
AVX (256-bit float), SSE4.2, POPCNT, F16C — but **no AVX2, no FMA, no BMI2**. In
`target-cpu` terms that is `x86-64-v2` + AVX = `-C target-cpu=ivybridge`.

- Set `RUSTFLAGS=-C target-cpu=ivybridge` for the two x86-64 desktop/Linux slices and the
  macOS Intel slice (via a per-target `[target.x86_64-*]` block in `.cargo/config.toml`,
  or per-job `RUSTFLAGS` in CI so the ARM/Apple-Silicon slices are unaffected).
- **Do not** set `x86-64-v3` — it emits AVX2/FMA that fault on the baseline.
- Any future hand-vectorized AVX2 path must use **runtime feature detection**
  (`is_x86_feature_detected!`) with an SSE/scalar fallback, never a compile-time `v3`
  floor. The binary's *baseline* stays Ivy Bridge; faster paths are dispatched at runtime.
- Pre-Ivy-Bridge CPUs (Sandy Bridge and older) will `SIGILL` — that is the **intended**
  floor, matching the stated "Ivy Bridge or newer" requirement.

---

## 4. Bundling the ancestry assets

Actual sizes in `~/.decodingus/ancestry/`:

| File | Size | Bundle? |
|------|------|---------|
| `ancestry_panel*` / `ancestry_pca*` (7 files) | ~9 MB total | **Yes** |
| `ancestry_freq_global` | 9.4 MB | **Yes** |
| `ancestry_manifest.json` | 4 KB | **Yes** |
| `genetic_map` | 39 MB | **Decision** |
| `ibd_panel` | 102 MB | **Decision** |

**Mechanism (recommended):** ship the asset files as **resources inside the image**
(macOS `.app/Contents/Resources/ancestry/`, Linux AppImage `usr/share`, Windows install
dir), and on first run **seed** any missing file into `~/.decodingus/ancestry/`, verifying
against the existing `AssetManifest` sha256 (the app already has `read_verified_asset` /
`ancestry_asset_status`). This keeps the runtime read path unchanged (still
`~/.decodingus/...`) and lets a later CDN download override a bundled asset transparently.

- Prefer **seeding data files** over `include_bytes!` — embedding 20 MB (let alone 160 MB)
  into the executable bloats every load and defeats incremental updates.
- The two big files: the app already has manifest-verified **CDN download** infra
  (`asset-manifest-verification`, panel pipeline → CDN publish). Two viable options:
  - **A — Full installer (CHOSEN):** bundle all 160 MB. Installer ≈ 180–200 MB. Everything
    works offline immediately; IBD ready out of the box. Simplest UX, heaviest download.
  - **B — Lean installer + lazy fetch:** bundle only the ~20 MB ancestry set; download
    `ibd_panel` + `genetic_map` (with manifest verification) on first IBD/segment use.
    Installer ≈ 40–60 MB. Lighter, but first IBD run needs network.

  **Decision (2026-06-16): A** — ship the full offline installer so Alpha testers have a
  zero-network, fully-working app including IBD. Revisit for GA if installer size becomes a
  friction point. The manifest-verified CDN download path stays in place as the update
  mechanism for refreshed assets.

---

## 5. Code signing & notarization (the operational crux)

| OS | Without signing | For Alpha | Cost |
|----|-----------------|-----------|------|
| macOS | Gatekeeper quarantine; user must right-click→Open, and notarization is required for a clean launch on current macOS | Get an **Apple Developer ID** ($99/yr), codesign + **notarize** the `.dmg`. cargo-packager/cargo-dist both automate this. | $99/yr |
| Windows | SmartScreen "unknown publisher" warning | Either ship unsigned with install instructions, or use **Azure Trusted Signing** (new, ~$10/mo, no EV dongle) or an OV cert | ~$10/mo+ |
| Linux | None needed | Ship AppImage + SHA256SUMS; optional detached GPG signature | Free |

**Decision (2026-06-16):** ship **macOS notarized** (Apple Developer ID, $99/yr — worth it
so testers don't fight Gatekeeper) and **Windows unsigned with a documented SmartScreen
work-around** for Alpha; add Windows signing (Azure Trusted Signing) before wider beta.
Linux: AppImage + `SHA256SUMS`, no signing.

---

## 6. Proposed CI shape (replaces `release.yml`)

Triggered on `v*` tags. One matrix job per image, plus a release-collection job.

- **macOS** (`macos-14`): add `x86_64-apple-darwin` target, build both slices with the
  Intel slice under `target-cpu=ivybridge`, `lipo` into a universal binary, `cargo-packager`
  → `.app`/`.dmg`, codesign + notarize with secrets.
- **Linux x86-64** (`ubuntu-22.04` or old-glibc container / zigbuild): `target-cpu=ivybridge`,
  `cargo-packager` → AppImage + `.deb` + `.rpm` + `.tar.gz`.
- **Linux ARM** (`ubuntu-24.04-arm`): native build, same images.
- **Windows x86-64** (`windows-latest`): `target-cpu=ivybridge`, `cargo-packager` → `.msi`
  + portable `.zip`; optional Authenticode.
- **Collect + release** (`softprops/action-gh-release`): attach all images + a
  `SHA256SUMS` file; mark prerelease for `*-alpha*`/`*-beta*`/`*-rc*` tags.

Asset seeding: a CI step stages the bundled ancestry files (fetched from the CDN by
manifest, not committed to git) into each image's resource dir before packaging.

---

## 7. Decisions & remaining open items

**Decided (2026-06-16):**
1. **Asset bundling:** Option A — full 160 MB offline installer.
2. **Signing:** macOS notarized (Apple Developer ID); Windows unsigned for Alpha, sign
   before beta; Linux checksums only.
3. **Packager:** `cargo-packager`.
4. **Implementation:** deferred — this doc stays research until the pipeline pass begins.

**Still open (defaults proposed, no blocker):**
5. **Windows ARM:** *defer past Alpha* unless an ARM-Windows tester appears.
6. **Linux glibc floor:** propose Debian 10 / glibc 2.28 via old-glibc container or
   `cargo-zigbuild`.

---

## References

- `docs/design/simd-optimization-targets.md` — the Ivy Bridge ISA constraint.
- `cargo-packager` — https://docs.rs/cargo-packager / CrabNebula.
- `cargo-dist` — https://axodotdev.github.io/cargo-dist/ (v0.31, maintained 2026).
- `cargo-zigbuild` — old-glibc / cross targeting via zig as linker.
