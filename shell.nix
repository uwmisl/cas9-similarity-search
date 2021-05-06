{
    nixpkgs ? import <nixpkgs> {}
,   lib ? nixpkgs.lib
,   python ? nixpkgs.python39
,   mach-nix
}:
let
    #requirements = ./requirements.txt;
    archiveQRMaker = nixpkgs.callPackage ./default.nix { pkgs=nixpkgs; inherit lib python mach-nix ; };

    photoshop_python_api = mach-nix.buildPythonPackage {
        src = builtins.fetchGit {
            url = "https://github.com/loonghao/photoshop-python-api";
            ref = "0.15.1";
        };
    };

    # Patch `packages.json` so that nix's *python* is used as default value for `python.pythonPath`.
    pathVsCodePythonHook = ''
        if [ -e ".vscode/settings.json" ]; then
            echo "Writing python.pythonPath..."
            sed 's|"python.pythonPath.*|"python.pythonPath": "${python}/yeet/bin/python"|g' -i ".vscode/settings.json"
        fi
    '';


in
nixpkgs.mkShell {
    #buildInputs=[
    #    nixpkgs.zsh
    #    nixpkgs.less
    #    nixpkgs.locale
    #    nixpkgs.ncurses
    #    python
    #];
    #propagatedBuildInputs=archiveQRMaker.bin.propagatedBuildInputs;

    propagatedBuildInputs=[ python ];

    shellHook=
       pathVsCodePythonHook +
    ''
       exec env zsh
    '' ;
}


