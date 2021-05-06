{
  description = "Code driving the UWMISL implementation of DNA similarity-search ";

  inputs = {
    nixpkgs = {
      url = "github:NixOS/nixpkgs";
    };

    flake-utils = {
      url = "github:numtide/flake-utils";
    };

    pypi-deps-db = {
      url = github:DavHau/pypi-deps-db;
      flake = false;
    };

    mach-nix = {
      url = github:DavHau/mach-nix;
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-utils.follows = "flake-utils";
      inputs.pypi-deps-db.follows = "pypi-deps-db";
    };
  };

  outputs = { self, nixpkgs, flake-utils, pypi-deps-db, mach-nix }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pname = "cas9-similarity-search";
        version = "0.5";
        pkgs = nixpkgs.legacyPackages.${system};

        mach-nix-utils = import mach-nix {
          inherit pkgs;
          pypiDataRev = pypi-deps-db.rev;
          pypiDataSha256 = pypi-deps-db.narHash;
        };

        python = mach-nix-utils.mkPython {
          python = "python27";
          requirements = builtins.readFile ./requirements.txt;
          packagesExtra = [
              pkgs.python
          ];
        };

        customOverrides = self: super: {
          # Overrides go here
        };

        app = mach-nix-utils.buildPythonApplication {
          inherit pname version;
          requirements = ''
            requests
            pillow
            tqdm
            h5py
            tables
            unireedsolomon
          '';
          # add missing dependencies whenever necessary.
          packagesExtra = [
            python
          ];
          src = ./.;
        };

        packageName = "cas9-similarity-search";
      in {
          packages.${packageName} = app;

          defaultPackage = self.packages.${system}.${packageName};

          devShell = import ./shell.nix { inherit python; nixpkgs=pkgs; mach-nix=mach-nix-utils; };
      });
}