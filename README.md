# dynatop

This R package coontains the core code to perform a dynamic TOPMODEL
simulation.

## TO DO
* Reformat example to reflect changes in debug example
* Provide R and Rcpp version of all functions
* clean out unused Rcpp functions
* Check function documentation and add missing variables
* What happens to cpp version of route_ex_eigen when eigen values of WV are
imaginary?
* check solution to initialisation when only one chanel
* change from channel to chanel
* Does channnel need vof??
* write checks
* add back in water balance?
* add back in output of states?

## Suggestions
* Move to object oriented?
* Check application of eigen vector routing?
* improve initialisation
* Change how initialisation is passed, currently loaded into workspace is
* allow checks to warn and change rather then fail
* Use RcppArmadillo - spped up matrix operation further + can use sparse matrices
* consider consistancy of sd and qsz
* qsz\_max initialised by atb\_bar rather then mean slope


## Have done

Core dynamic TOPMODEL code
* Gone through Peter code to find functions actually beign used -
  stored as 'clean' branch
* Altered input files to allow for more easy reparameterisation and
  use of multiple inputs (used named columns, parameters not implicit
  by columns order)
* rewritten initialisation and main execution loop of dynamic TOPMODEL
  so that
  * Chanels are handled explicily not as 'special cases'
    (parameterisations) of hillslope units
  * Outputs can be selected

## Hydrological Response Units

Three Hydrological Response Units (HRU) ypes are defined.

### Hillslope

States [all m]

* Excess water on surface - sxs
* Root zone stroage - srz
* unsaturated zone storage - suz
* storage deficit in saturated zone - sd

Inputs

* precipitation - precip [m/hr]
* potential evapotranspiration - pet [m/hr]

Outputs

* base flow discharge
* Surface excess discharge
* 

Internal Fluxes


Tempory states

qin qbf quz qof
* Flow from root zone to 


### Channel_Output

States

None

Inputs

* Fluxes

Conditions


### Outlet

RU 


**Edit a file, create a new file, and clone from Bitbucket in under 2 minutes**

When you're done, you can delete the content in this README and update the file with details for others getting started with your repository.

*We recommend that you open this README in another tab as you perform the tasks below. You can [watch our video](https://youtu.be/0ocf7u76WSo) for a full demo of all the steps in this tutorial. Open the video in a new tab to avoid leaving Bitbucket.*

---

## Edit a file

You’ll start by editing this README file to learn how to edit a file in Bitbucket.

1. Click **Source** on the left side.
2. Click the README.md link from the list of files.
3. Click the **Edit** button.
4. Delete the following text: *Delete this line to make a change to the README from Bitbucket.*
5. After making your change, click **Commit** and then **Commit** again in the dialog. The commit page will open and you’ll see the change you just made.
6. Go back to the **Source** page.

---

## Create a file

Next, you’ll add a new file to this repository.

1. Click the **New file** button at the top of the **Source** page.
2. Give the file a filename of **contributors.txt**.
3. Enter your name in the empty file space.
4. Click **Commit** and then **Commit** again in the dialog.
5. Go back to the **Source** page.

Before you move on, go ahead and explore the repository. You've already seen the **Source** page, but check out the **Commits**, **Branches**, and **Settings** pages.

---

## Clone a repository

Use these steps to clone from SourceTree, our client for using the repository command-line free. Cloning allows you to work on your files locally. If you don't yet have SourceTree, [download and install first](https://www.sourcetreeapp.com/). If you prefer to clone from the command line, see [Clone a repository](https://confluence.atlassian.com/x/4whODQ).

1. You’ll see the clone button under the **Source** heading. Click that button.
2. Now click **Check out in SourceTree**. You may need to create a SourceTree account or log in.
3. When you see the **Clone New** dialog in SourceTree, update the destination path and name if you’d like to and then click **Clone**.
4. Open the directory you just created to see your repository’s files.

Now that you're more familiar with your Bitbucket repository, go ahead and add a new file locally. You can [push your change back to Bitbucket with SourceTree](https://confluence.atlassian.com/x/iqyBMg), or you can [add, commit,](https://confluence.atlassian.com/x/8QhODQ) and [push from the command line](https://confluence.atlassian.com/x/NQ0zDQ).
