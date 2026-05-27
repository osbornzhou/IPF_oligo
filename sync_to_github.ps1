param(
    [string]$Owner = "osbornzhou",
    [string]$Repo = "IPF_oligo",
    [string]$Branch = "main",
    [string]$Token = $env:GITHUB_TOKEN
)

if (-not $Token) {
    throw "Set GITHUB_TOKEN to a fine-grained GitHub token with Contents: Read and write permission, then rerun this script."
}

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Files = @(
    "README.md",
    "scripts/README.md",
    "manuscript/repository_deposit_README.md",
    "manuscript/ipf_oligo_ml_bmc_genomics_manuscript_draft.md",
    "manuscript/ipf_oligo_ml_bmc_genomics_submission_ready.docx",
    "manuscript/ipf_oligo_ml_bmc_genomics_submission_ready.pdf"
)

$Files += Get-ChildItem -Path (Join-Path $Root "manuscript/additional_files") -File | ForEach-Object {
    "manuscript/additional_files/$($_.Name)"
}
$Files += Get-ChildItem -Path (Join-Path $Root "manuscript/reproducibility_package") -File | ForEach-Object {
    "manuscript/reproducibility_package/$($_.Name)"
}

$Headers = @{
    Authorization = "Bearer $Token"
    Accept = "application/vnd.github+json"
    "X-GitHub-Api-Version" = "2022-11-28"
}

foreach ($Path in $Files) {
    $LocalPath = Join-Path $Root ($Path -replace "/", [System.IO.Path]::DirectorySeparatorChar)
    if (-not (Test-Path -LiteralPath $LocalPath)) {
        Write-Warning "Missing local file: $Path"
        continue
    }

    $UriPath = [System.Uri]::EscapeDataString($Path).Replace("%2F", "/")
    $Uri = "https://api.github.com/repos/$Owner/$Repo/contents/$UriPath"
    $Sha = $null
    try {
        $Existing = Invoke-RestMethod -Method Get -Uri "$Uri`?ref=$Branch" -Headers $Headers
        $Sha = $Existing.sha
    } catch {
        $Sha = $null
    }

    $Content = [Convert]::ToBase64String([System.IO.File]::ReadAllBytes($LocalPath))
    $Body = @{
        message = "Synchronize submission package reproducibility files"
        content = $Content
        branch = $Branch
    }
    if ($Sha) {
        $Body.sha = $Sha
    }

    $JsonBody = $Body | ConvertTo-Json -Depth 5
    Invoke-RestMethod -Method Put -Uri $Uri -Headers $Headers -Body $JsonBody -ContentType "application/json" | Out-Null
    Write-Host "Uploaded $Path"
}

Write-Host "GitHub synchronization completed for $Owner/$Repo on branch $Branch."
