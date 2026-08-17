$ErrorActionPreference = 'Stop'
try {
    $Root = $PSScriptRoot
    $AppRoot = Join-Path $Root 'app'
    if (-not (Test-Path -LiteralPath (Join-Path $AppRoot 'src'))) {
        throw "Application files not found: $AppRoot"
    }

    $PortablePythonw = Join-Path $Root 'python\pythonw.exe'
    $PortablePython = Join-Path $Root 'python\python.exe'
    if (Test-Path -LiteralPath $PortablePythonw -PathType Leaf) {
        $pythonw = $PortablePythonw
    } elseif (Test-Path -LiteralPath $PortablePython -PathType Leaf) {
        $pythonw = $PortablePython
    } else {
        throw "Portable Python not found. Keep the bundled 'python' folder beside this launcher."
    }

    $env:PYTHONPATH = Join-Path $AppRoot 'src'
    $PortableIgBlastBin = Join-Path $Root 'tools\igblast-1.21.0\bin'
    $PortableIgBlast = Join-Path $PortableIgBlastBin 'igblastn.exe'
    if (-not (Test-Path -LiteralPath $PortableIgBlast -PathType Leaf)) {
        throw "Portable IgBLAST not found: $PortableIgBlast"
    }
    $env:PATH = "$PortableIgBlastBin;$env:PATH"

    $RefdataRoot = Join-Path $Root 'refdata\IgBlast_refdata_edit_imgt'
    foreach ($DbName in @('IMGT_IGHV.imgt', 'IMGT_IGHD.imgt', 'IMGT_IGHJ.imgt')) {
        $Prefix = Join-Path (Join-Path $RefdataRoot 'db') $DbName
        if (-not (
            (Test-Path -LiteralPath "$Prefix.nsq" -PathType Leaf) -or
            (Test-Path -LiteralPath "$Prefix.njs" -PathType Leaf) -or
            (Test-Path -LiteralPath "$Prefix.nin" -PathType Leaf)
        )) {
            throw "IgBLAST database components not found for prefix: $Prefix"
        }
    }
    $Auxiliary = Join-Path $RefdataRoot 'optional_file\human_gl.aux'
    if (-not (Test-Path -LiteralPath $Auxiliary -PathType Leaf)) {
        throw "IgBLAST auxiliary file not found: $Auxiliary"
    }

    Start-Process -FilePath $pythonw -ArgumentList @('-m', 'airr_igblast_paired', 'gui') -WorkingDirectory $AppRoot
} catch {
    Add-Type -AssemblyName PresentationFramework
    [System.Windows.MessageBox]::Show($_.Exception.Message, 'CPM Paired Fastq IgBLAST AIRR tsv Ver3.0') | Out-Null
}
