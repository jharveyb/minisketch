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
    in
    {
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
          coz
          python3

          # Debugging / profiling / benchmarking (tools/bench_ab.sh et al.)
          gdb
          valgrind
          hyperfine
          samply
          perf
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
