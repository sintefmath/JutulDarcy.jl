module.exports = async ({ github, context, core }) => {
  const { owner, repo } = context.repo;
  const pr = context.payload.pull_request;
  // A rerun for an old commit must not claim that it documents the current head.
  const { data: current } = await github.rest.pulls.get({ owner, repo, pull_number: pr.number });
  if (current.state !== 'open' || current.head.sha !== pr.head.sha) return;
  const url = `https://${owner}.github.io/${repo}/previews/PR${pr.number}/`;
  const marker = '<!-- jutuldarcy-doc-preview -->';
  const body = `${marker}\nDocumentation preview: [Open preview](${url})\n\nBuilt from ${pr.head.sha}.`;
  const comments = await github.paginate(github.rest.issues.listComments, {
    owner, repo, issue_number: pr.number, per_page: 100,
  });
  const previous = comments.find(c => c.user?.login === 'github-actions[bot]' && c.body?.startsWith(marker));
  if (previous) {
    await github.rest.issues.updateComment({ owner, repo, comment_id: previous.id, body });
  } else {
    await github.rest.issues.createComment({ owner, repo, issue_number: pr.number, body });
  }
  await core.summary.addLink('Documentation preview', url).write();
};
