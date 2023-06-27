# Mangement of the GitHub project.

resource "github_repository" "varfish-server-worker" {
  name        = "varfish-server-worker"
  description = "Rust-based background worker for varfish-server"
  visibility  = "public"

  has_issues      = true
  has_wiki        = true
  has_downloads   = true
  has_discussions = false
  has_projects    = false

  allow_auto_merge   = true
  allow_rebase_merge = false
  allow_merge_commit = false

  delete_branch_on_merge = true

  vulnerability_alerts = true

  squash_merge_commit_message = "BLANK"
  squash_merge_commit_title   = "PR_TITLE"
  merge_commit_message        = "PR_BODY"
  merge_commit_title          = "PR_TITLE"
}
