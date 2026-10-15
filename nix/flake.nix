{
  description = "coz and mull packaged nixpkgs-style for the minisketch devshell";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-26.05";

  outputs =
    { self, nixpkgs }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
      ];
      forAllSystems =
        f:
        nixpkgs.lib.genAttrs systems (
          system:
          f (
            import nixpkgs {
              inherit system;
              overlays = [ self.overlays.default ];
            }
          )
        );
    in
    {
      overlays.default = final: prev: {
        coz = final.callPackage ./pkgs/by-name/co/coz/package.nix { };
        mull = final.callPackage ./pkgs/by-name/mu/mull/package.nix { };
      };

      packages = forAllSystems (pkgs: {
        inherit (pkgs) coz mull;
      });

      devShells = forAllSystems (pkgs: {
        default = pkgs.mkShell {
          packages = [
            pkgs.coz
            pkgs.mull
            pkgs.mull.llvmPackages.clang
          ];
          MULL_PLUGIN = pkgs.mull.plugin;
        };
      });

      formatter = forAllSystems (pkgs: pkgs.nixfmt-tree);
    };
}
