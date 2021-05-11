{
    nixpkgs ? import <nixpkgs> {}
,   lib ? nixpkgs.lib
,   python ? nixpkgs.python39
,   mach-nix
}:
let
    # Patch `packages.json` so that nix's *python* is used as default value for `python.pythonPath`.
    pathVsCodePythonHook = ''
        if [ -e ".vscode/settings.json" ]; then
            echo "Writing python.pythonPath..."
            sed 's|"python.pythonPath.*|"python.pythonPath": "${python}/yeet/bin/python"|g' -i ".vscode/settings.json"
        fi
    '';

in
nixpkgs.mkShell {
    propagatedBuildInputs=[ python ];

    shellHook=
       #pathVsCodePythonHook +
    ''
       exec env zsh
    '' ;
}


