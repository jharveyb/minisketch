{
  lib,
  stdenv,
  fetchurl,
  autoPatchelfHook,
  dpkg,
  installShellFiles,
  llvmPackages,
}:

let
  llvmMajor = lib.versions.major llvmPackages.llvm.version;

  # Upstream builds with Bazel using hermetic toolchains, so we repackage the
  # official per-LLVM-version .deb releases instead.
  sources = {
    "21-x86_64-linux" = {
      llvmVersion = "21.1.8";
      arch = "amd64";
      hash = "sha256-HwrpIcmRqqf3vEA31G5hLpEL91/tcisJABK1qH6bK04=";
    };
    "21-aarch64-linux" = {
      llvmVersion = "21.1.8";
      arch = "aarch64";
      hash = "sha256-Snl1JxN54idLImN0IjTuPrrhJixpQ1Lr/nj6dlPEHss=";
    };
  };

  source =
    sources."${llvmMajor}-${stdenv.hostPlatform.system}"
      or (throw "mull: no release for LLVM ${llvmMajor} on ${stdenv.hostPlatform.system}");
in
stdenv.mkDerivation (finalAttrs: {
  pname = "mull";
  version = "0.34.1";

  src = fetchurl {
    url = "https://github.com/mull-project/mull/releases/download/${finalAttrs.version}/Mull-${llvmMajor}-${finalAttrs.version}-LLVM-${source.llvmVersion}-ubuntu-${source.arch}-26.04.deb";
    inherit (source) hash;
  };

  nativeBuildInputs = [
    autoPatchelfHook
    dpkg
    installShellFiles
  ];

  buildInputs = [
    llvmPackages.libllvm.lib
    llvmPackages.libclang.lib
    (lib.getLib stdenv.cc.cc)
  ];

  installPhase = ''
    runHook preInstall

    install -Dm755 -t $out/bin usr/bin/*
    install -Dm755 -t $out/lib usr/lib/mull-ir-frontend-${llvmMajor}
    installManPage usr/share/man/man1/*

    for tool in runner reporter; do
      installShellCompletion \
        --cmd mull-$tool-${llvmMajor} \
        --bash usr/share/bash-completion/completions/mull-$tool-${llvmMajor} \
        --fish usr/share/fish/vendor_completions.d/mull-$tool-${llvmMajor}.fish \
        --zsh usr/share/zsh/vendor-completions/_mull-$tool-${llvmMajor}
    done

    for tool in instrument runner reporter; do
      ln -s mull-$tool-${llvmMajor} $out/bin/mull-$tool
    done

    runHook postInstall
  '';

  passthru = {
    inherit llvmPackages;
    # Use with a clang of the same major version:
    #   clang -fpass-plugin=${mull.plugin} -g -grecord-command-line
    plugin = "${finalAttrs.finalPackage}/lib/mull-ir-frontend-${llvmMajor}";
  };

  meta = {
    description = "Practical mutation testing and fault injection for C and C++";
    homepage = "https://github.com/mull-project/mull";
    changelog = "https://github.com/mull-project/mull/releases/tag/${finalAttrs.version}";
    license = lib.licenses.asl20;
    sourceProvenance = [ lib.sourceTypes.binaryNativeCode ];
    mainProgram = "mull-runner-${llvmMajor}";
    platforms = [
      "x86_64-linux"
      "aarch64-linux"
    ];
    maintainers = [ ];
  };
})
