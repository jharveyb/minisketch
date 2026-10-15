{
  lib,
  stdenv,
  fetchFromGitHub,
  cmake,
  lief,
  python3,
}:

let
  # Upstream's CMake fetches this fork (branch improved_dwarf5_support) for
  # DWARF5 support, which nixpkgs' libelfin (aclements upstream) lacks.
  libelfin-src = fetchFromGitHub {
    owner = "plasma-umass";
    repo = "libelfin";
    rev = "5a032d5896f98aa57200e6e58e93f1dac0527a90";
    hash = "sha256-nIv2z5f8IEUfPxstCHm0tGviM+Nz3zKB8S2GribWdLY=";
  };
in
stdenv.mkDerivation (finalAttrs: {
  pname = "coz";
  version = "0.2.5";

  src = fetchFromGitHub {
    owner = "plasma-umass";
    repo = "coz";
    tag = "v${finalAttrs.version}";
    hash = "sha256-h28z3Tt4CIHKcxyzwYM3wL7D3/4EvJ1q3PCfkEfcFi4=";
  };

  postPatch = ''
    # Neither package metadata nor pyproject.toml is available at runtime
    substituteInPlace coz \
      --replace-fail "__version__ = _read_version()" "__version__ = '${finalAttrs.version}'"
  '';

  nativeBuildInputs = [ cmake ];

  buildInputs = [
    lief
    python3
  ];

  cmakeFlags = [
    (lib.cmakeBool "USE_SYSTEM_LIEF" true)
    (lib.cmakeFeature "FETCHCONTENT_SOURCE_DIR_LIBELFIN" "${libelfin-src}")
    (lib.cmakeBool "FETCHCONTENT_FULLY_DISCONNECTED" true)
  ];

  # Tests need perf_event_open, which is unavailable in the build sandbox
  doCheck = false;

  meta = {
    homepage = "https://github.com/plasma-umass/coz";
    changelog = "https://github.com/plasma-umass/coz/releases/tag/v${finalAttrs.version}";
    description = "Profiler based on casual profiling";
    mainProgram = "coz";
    license = lib.licenses.bsd2;
    platforms = lib.platforms.linux;
    maintainers = with lib.maintainers; [
      zimbatm
      aleksana
    ];
  };
})
