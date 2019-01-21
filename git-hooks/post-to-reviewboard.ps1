param([string]$range)
If ($range.Contains('..')) {
    echo "Submitting to RBCommon using RBTool git range:"
    echo $range
} Else {
    echo 'This does not appear to be range exiting'
    exit 1
}

$commits = git rev-list $range --topo-order --reverse
foreach ($newrev in $commits) {
    echo "Submitting revision:"
    echo $newrev
    #get first name from github commit email address
    $githubname=$(git show --format="%ae" $newrev --no-patch)
    echo "Looking up RB username"
    $rbUserList = @{}
    $rbUserList.add('lkopp@akoyabio.com', 'lorenzo')
    $rbUserList.add('adubey@akoyabio.com', 'ajay')
    $rbUserList.add('aostapenko@akoyabio.com', 'alex')
    $rbUserList.add('ccoltharp@akoyabio.com', 'carla')
    $rbUserList.add('emarquez@akoyabio.com', 'elbano')
    $rbUserList.add('eedwards@akoyabio.com', 'elizabeth')
    $rbUserList.add('kjohnson@akoyabio.com', 'kent')
    $rbUserList.add('pmiller@akoyabio.com', 'peter')
    $rbUserList.add('ssheedy@akoyabio.com', 'stephen')
    $rbUserList.add('wzhang@akoyabio.com', 'wen')
    #add the legacy emails for old commits
    $rbUserList.add('lorenzo.kopp@perkinelmer.com', 'lorenzo')
    $rbUserList.add('ajay.dubey@perkinelmer.com', 'ajay')
    $rbUserList.add('alex.ostapenko@perkinelmer.com', 'alex')
    $rbUserList.add('elbano.marquez@perkinelmer.com', 'elbano')
    $rbUserList.add('elizabeth.edwards@perkinelmer.com', 'elizabeth')
    $rbUserList.add('kent.johnson@perkinelmer.com', 'kent')
    $rbUserList.add('peter.miller@perkinelmer.com', 'peter')
    $rbUserList.add('stephen.sheedy@perkinelmer.com', 'stephen')

    $first = $rbUserList[$githubname]
    echo "  First name found: $first"

    #get ticket number from commit message
    $gitsubject=$(git show --format="%s" $newrev --no-patch)
    $found=$gitsubject -match '(?<=ticket #)\d+'
    If ($found) {
      $bugid=$matches[0]
      $bugstring="--bugs-closed=$bugid"
      echo "  Ticket id found: $bugid"
    } Else {
      $bugstring=""
      echo "  No ticket id found."
    }

    #get branch name
    $branchname=$(git rev-parse --abbrev-ref HEAD)

    #get random target to review post
    DO {
      $targetpeople=Get-Random -Maximum ('lorenzo', 'stephen', 'elizabeth', 'peter', 'wen', 'kent')
    } While ($targetpeople -eq $first)
    echo "  Random target reviewer: $targetpeople"

    #setup command and post it NOTE: add a -d for debug mode if things break!
    $reviewboardresponse = $(rbt post -p --submit-as=$first --username=maxwell --password=maxwell $bugstring --branch $branchname --target-people=$targetpeople $newrev) | Select-String -Pattern 'Review request'
    $gitshortid = $newrev.SubString(0,8)
    echo "  RB Response:"$reviewboardresponse
    $found = $reviewboardresponse -match 'Review request #(\d+)'
    $reviewboardid = $matches[1]
    If ($found) {
      echo "  GitHub:$gitshortid (reviewboard:$reviewboardid)"
    } Else {
      echo "  GitHub:$gitshortid (reviewboard:NOTFOUND)"
    }
}
exit 0

#EXAMPLE: rbt post -p -d --submit-as=lorenzo --username=maxwell --password=maxwell --bugs-closed=2 --branch master --target-people=lorenzo e28e7f12fa34a7f598961b28f5e6b3d08f92b48f

# -d debug
# -p - publish the draft
# --submit-as USERNAME - to set the sender
# --username - set maxwell to update RB
# --password - set maxwell to update RB
# --bugs-closed BUG_ID - set to the ticket ids
# --branch - freeform only doesn't effect diff
# --target-people USERNAME - the person doing the review should be random - self
# --parent BRANCH - the real branch name to diff from
# --revision-range REV1[:REV2] - sets a diff by actual revision <- we want this one
