{
  description = "minisketch dev environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-26.05";
    # coz and mull, packaged nixpkgs-style (pkgs/by-name + overlay). Follows
    # our nixpkgs: its mull/coz hashes were verified against nixos-26.05.
    minisketch-deps = {
      url = "path:./nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  # x86_64-linux only: the dev VM (CLMUL/AVX2 passthrough) is the one target.
  # nix/ also builds for aarch64-linux if that is ever needed here.
  outputs = { self, nixpkgs, minisketch-deps }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs {
        inherit system;
        overlays = [ minisketch-deps.overlays.default ];
      };

      # The one clang in the shell: mull's IR plugin is LLVM-version-matched,
      # and mull.llvmPackages is the set its .deb was built against.
      mullClang = pkgs.mull.llvmPackages.clang;

      # MCP server behind .mcp.json's `bkb` entry (Bitcoin Knowledge Base,
      # the skill in .claude/skills/bitcoin-knowledge-base). Same pin as
      # ~/gossip_observer/flake.nix; .mcp.json runs it via `nix run .#bkb-mcp`
      # so it works without being on the PATH Claude Code was started with.
      bkb-mcp = pkgs.rustPlatform.buildRustPackage {
        pname = "bkb-mcp";
        version = "0.1.0";
        src = pkgs.fetchFromGitHub {
          owner = "tnull";
          repo = "bitcoin-knowledge-base";
          rev = "3d917d6aa61d8ea5bfab0f3dadfc9a1173bd8b85";
          hash = "sha256-rRwo+wWksq7mTcbo6EoTzcb3+vE8wU7Bj0cmiWO95I0=";
        };
        cargoHash = "sha256-4kUq0LZjrVQrUIzAmvJpxNJkxjeP6puvSZL1xZU0tPc=";
        cargoBuildFlags = [ "-p" "bkb-mcp" ];
        doCheck = false;
        nativeBuildInputs = [ pkgs.pkg-config ];
        buildInputs = [ pkgs.openssl ];
        meta.mainProgram = "bkb-mcp";
      };
    in
    {
      packages.${system} = { inherit bkb-mcp; };

      devShells.${system}.default = pkgs.mkShell {
        packages = with pkgs; [
          # Build chain (autotools + CMake dual build). gcc is listed before
          # mullClang so the generic `cc`/`c++` names resolve to gcc; clang is
          # only reached by name (clang, clang++).
          autoconf
          automake
          libtool
          pkg-config
          gcc
          gnumake
          cmake

          # Mutation testing / causal profiling (overlay). tools/mutation.sh
          # picks up $MULL_PLUGIN and clang++ from here; tools/coz.sh needs
          # coz.h (CMake MINISKETCH_COZ) and python3 for coz_summary.py.
          mull
          mullClang
          # clangd for the clangd-lsp Claude plugin / editors: same LLVM set,
          # not a second toolchain. Reads build/compile_commands.json.
          pkgs.mull.llvmPackages.clang-tools
          coz
          python3

          # Debugging / profiling / benchmarking (tools/bench_ab.sh et al.)
          gdb
          valgrind
          hyperfine
          samply
          perf

          # scripts/claude_setup.sh
          jq
        ];

        # Libraries: buildInputs get -isystem include dirs / CMake prefix paths.
        #  - boost: Boost.Test for the CMake unit-tests target (also the Mull
        #    killer binary); src/CMakeLists.txt skips it with only a warning.
        #  - coz: coz.h for -DMINISKETCH_COZ=ON.
        buildInputs = with pkgs; [
          boost
          coz
        ];

        env = {
          # Usage: clang++ -fpass-plugin=$MULL_PLUGIN -g -grecord-command-line
          MULL_PLUGIN = "${pkgs.mull.plugin}";
          # CMake's default for new build dirs; clangd finds
          # build/compile_commands.json on its own.
          CMAKE_EXPORT_COMPILE_COMMANDS = "ON";
          # compile_commands.json names the nix g++ wrapper, which clangd will
          # not run to learn the libstdc++ include paths unless allowed (the
          # clangd-lsp plugin's args are fixed, so pass it via the env).
          CLANGD_FLAGS = "--query-driver=/nix/store/*-gcc-wrapper-*/bin/g++";
        };

        # mull.yml is located through MULL_CONFIG (CMake compiles from the
        # build dir, so mull's cwd lookup misses it). tools/mutation.sh sets it
        # itself; exporting it here covers ad-hoc instrumented builds too.
        #
        # mullClang's setup hook also sets CC/CXX to clang, last hook wins;
        # keep gcc the default compiler (as CI's default is g++), clang opt-in.
        shellHook = ''
          export CC=gcc CXX=g++
          if [ -z "''${MULL_CONFIG:-}" ] && [ -f "$PWD/mull.yml" ]; then
            export MULL_CONFIG="$PWD/mull.yml"
          fi
        '';
      };
    };
}
